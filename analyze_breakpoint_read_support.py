#!/usr/bin/env python3
"""Independent long-read validation of two qHTSF4.1 breakpoints."""
import argparse, csv, os, re, sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
try:
    import pysam
except ImportError:
    sys.exit("ERROR: pysam is required: conda install -c bioconda pysam")

@dataclass
class Dataset:
    reference:str; sample:str; bam:str; fasta:str; chrom:str; start:int; end:int

def load_config(path):
    need=["reference","sample","bam","fasta","chrom","region_start","region_end"]
    p=Path(path).resolve()
    if not p.is_file(): raise RuntimeError(f"Config not found: {p}")
    out=[]
    with p.open(newline="") as f:
        r=csv.DictReader(f,delimiter="\t")
        missing=[x for x in need if not r.fieldnames or x not in r.fieldnames]
        if missing: raise RuntimeError("Config missing columns: "+", ".join(missing))
        for n,row in enumerate(r,2):
            if not any((v or "").strip() for v in row.values()): continue
            try: start,end=int(row["region_start"]),int(row["region_end"])
            except ValueError: raise RuntimeError(f"Config line {n}: coordinates must be integers")
            if start>=end: raise RuntimeError(f"Config line {n}: region_start must be < region_end")
            bam,fa=Path(row["bam"]),Path(row["fasta"])
            if not bam.is_absolute(): bam=p.parent/bam
            if not fa.is_absolute(): fa=p.parent/fa
            out.append(Dataset(row["reference"],row["sample"],str(bam.resolve()),str(fa.resolve()),row["chrom"],start,end))
    if not out: raise RuntimeError("Config contains no datasets")
    return out

def validate(ds):
    for label,path in (("BAM",ds.bam),("FASTA",ds.fasta)):
        if not os.path.isfile(path) or os.path.getsize(path)==0: raise RuntimeError(f"{ds.reference}: missing/empty {label}: {path}")
    pysam.quickcheck(ds.bam); bam=pysam.AlignmentFile(ds.bam,"rb"); fa=pysam.FastaFile(ds.fasta)
    if not bam.has_index(): raise RuntimeError(f"{ds.reference}: BAM index missing")
    if ds.chrom not in bam.references or ds.chrom not in fa.references: raise RuntimeError(f"{ds.reference}: {ds.chrom} absent from BAM/FASTA")
    if bam.get_reference_length(ds.chrom)!=fa.get_reference_length(ds.chrom): raise RuntimeError(f"{ds.reference}: BAM/FASTA chromosome length mismatch")
    length=fa.get_reference_length(ds.chrom)
    if not 1<=ds.start<ds.end<=length: raise RuntimeError(f"{ds.reference}: interval outside chromosome")
    return bam,fa,length

def strand(r): return "-" if r.is_reverse else "+"
def atype(r): return "secondary" if r.is_secondary else ("supplementary" if r.is_supplementary else "primary")

def clips(r):
    if not r.cigartuples: return 0,0
    left=r.cigartuples[0][1] if r.cigartuples[0][0] in (4,5) else 0
    right=r.cigartuples[-1][1] if r.cigartuples[-1][0] in (4,5) else 0
    return left,right

def qinterval(r):
    qlen=r.infer_read_length() or r.query_length or 0; a,b=r.query_alignment_start,r.query_alignment_end
    hl=hr=0
    if r.cigartuples:
        hl=r.cigartuples[0][1] if r.cigartuples[0][0]==5 else 0
        hr=r.cigartuples[-1][1] if r.cigartuples[-1][0]==5 else 0
    a+=hl; b+=hl; qlen+=hl+hr
    if r.is_reverse and qlen: a,b=qlen-b,qlen-a
    return a+1,b

def cigar_ref_len(cigar):
    return sum(int(n) for n,op in re.findall(r"(\d+)([MIDNSHP=X])",cigar) if op in "MDN=X")

def parse_sa(r):
    if not r.has_tag("SA"): return []
    out=[]
    for item in r.get_tag("SA").rstrip(";").split(";"):
        f=item.split(",")
        if len(f)!=6: continue
        try:
            start=int(f[1]); out.append(dict(chrom=f[0],start=start,end=start+cigar_ref_len(f[3])-1,strand=f[2],cigar=f[3],mapq=int(f[4]),nm=int(f[5])))
        except ValueError: pass
    return out

def find_n_runs(fa,chrom,lo,hi,min_len):
    seq=fa.fetch(chrom,lo-1,hi).upper()
    return [(lo+m.start(),lo+m.end()-1,m.end()-m.start()) for m in re.finditer(r"N{%d,}"%min_len,seq)]

def distance(point,start,end):
    return 0 if start<=point<=end else (start-point if point<start else point-end)

def cluster(events,tol):
    groups=[]
    for e in sorted(events,key=lambda x:x["coordinate"]):
        if not groups or e["coordinate"]-groups[-1]["max"]>tol:
            groups.append(dict(events=[e],min=e["coordinate"],max=e["coordinate"]))
        else:
            groups[-1]["events"].append(e); groups[-1]["max"]=e["coordinate"]
    out=[]
    for g in groups:
        ids=sorted({e["read_id"] for e in g["events"]})
        coords=sorted({e["coordinate"] for e in g["events"]})
        peak=sorted(coords,key=lambda c:(-len({e["read_id"] for e in g["events"] if e["coordinate"]==c}),c))[0]
        out.append(dict(cluster_start=g["min"],cluster_end=g["max"],peak_coordinate=peak,unique_reads=len(ids),event_count=len(g["events"]),read_ids=",".join(ids)))
    return sorted(out,key=lambda x:(-x["unique_reads"],x["cluster_start"]))

def analyse(ds,a):
    bam,fa,chrom_len=validate(ds); boundaries=[("LEFT",ds.start),("RIGHT",ds.end)]
    local=defaultdict(list); names=set()
    for label,bp in boundaries:
        lo,hi=max(1,bp-a.window),min(chrom_len,bp+a.window)
        for r in bam.fetch(ds.chrom,lo-1,hi):
            if r.is_unmapped or r.is_secondary: continue
            local[label].append(r); names.add(r.query_name)
    allrecs=defaultdict(list); bam.reset()
    for r in bam.fetch(until_eof=True):
        if r.query_name in names and not r.is_unmapped: allrecs[r.query_name].append(r)
    summary=[]; clipped=[]; sa_rows=[]; endpoint=[]; gaps=[]; alignments=[]; clipped_names=set()
    for label,bp in boundaries:
        good=[r for r in local[label] if r.reference_name==ds.chrom and r.mapping_quality>=a.min_mapq and r.query_alignment_length>=a.min_anchor]
        spanning={r.query_name for r in good if r.reference_start+1<=bp-a.span_anchor and r.reference_end>=bp+a.span_anchor}
        events=[]; clip_ids=set(); sa_ids=set()
        for r in good:
            lc,rc=clips(r)
            for side,clip_bp,coord in (("REF_LEFT",lc,r.reference_start+1),("REF_RIGHT",rc,r.reference_end)):
                if clip_bp<a.min_clip or abs(coord-bp)>a.window: continue
                clip_ids.add(r.query_name); clipped_names.add(r.query_name); events.append(dict(coordinate=coord,read_id=r.query_name))
                sas=parse_sa(r)
                if sas: sa_ids.add(r.query_name)
                qa,qb=qinterval(r)
                clipped.append(dict(reference=ds.reference,breakpoint_label=label,configured_breakpoint=bp,read_id=r.query_name,
                    alignment_type=atype(r),alignment_location=f"{ds.chrom}:{r.reference_start+1}-{r.reference_end}",alignment_strand=strand(r),
                    endpoint_side=side,endpoint_coordinate=coord,offset_from_configured_bp=coord-bp,clip_bp=clip_bp,mapq=r.mapping_quality,
                    cigar=r.cigarstring,query_interval=f"{qa}-{qb}",sa_entry_count=len(sas)))
                for sa in sas:
                    near=[]
                    if sa["chrom"]==ds.chrom:
                        if min(abs(sa["start"]-ds.start),abs(sa["end"]-ds.start))<=a.window: near.append("LEFT")
                        if min(abs(sa["start"]-ds.end),abs(sa["end"]-ds.end))<=a.window: near.append("RIGHT")
                    sa_rows.append(dict(reference=ds.reference,breakpoint_label=label,configured_breakpoint=bp,read_id=r.query_name,
                        source_location=f"{ds.chrom}:{r.reference_start+1}-{r.reference_end}",source_type=atype(r),source_strand=strand(r),source_cigar=r.cigarstring,
                        sa_chrom=sa["chrom"],sa_start=sa["start"],sa_end=sa["end"],sa_strand=sa["strand"],sa_cigar=sa["cigar"],
                        sa_mapq=sa["mapq"],sa_nm=sa["nm"],sa_near_boundary="+".join(near) if near else "NONE"))
        clusters=sorted(cluster(events,a.cluster_tolerance),
                        key=lambda x:(-x["unique_reads"],abs(x["peak_coordinate"]-bp),x["cluster_start"]))
        for rank,c in enumerate(clusters,1):
            endpoint.append(dict(reference=ds.reference,breakpoint_label=label,configured_breakpoint=bp,cluster_rank=rank,**c,peak_offset=c["peak_coordinate"]-bp))
        best=clusters[0] if clusters else None
        lo,hi=max(1,bp-a.window),min(chrom_len,bp+a.window); nruns=find_n_runs(fa,ds.chrom,lo,hi,a.min_n_run)
        nearest=min(nruns,key=lambda x:distance(bp,x[0],x[1])) if nruns else None; gap_union=set()
        for gs,ge,glen in nruns:
            span={r.query_name for r in good if r.reference_start+1<=gs-a.span_anchor and r.reference_end>=ge+a.span_anchor}
            gap_union|=span
            gaps.append(dict(reference=ds.reference,breakpoint_label=label,configured_breakpoint=bp,gap_start=gs,gap_end=ge,gap_length=glen,
                distance_from_breakpoint=distance(bp,gs,ge),breakpoint_inside_gap="YES" if gs<=bp<=ge else "NO",
                gap_spanning_read_count=len(span),gap_spanning_read_ids=",".join(sorted(span))))
        summary.append(dict(reference=ds.reference,read_sample=ds.sample,chrom=ds.chrom,breakpoint_label=label,configured_breakpoint=bp,
            continuous_spanning_reads=len(spanning),continuous_read_ids=",".join(sorted(spanning)),clipped_reads=len(clip_ids),clipped_read_ids=",".join(sorted(clip_ids)),
            clipped_reads_with_sa=len(sa_ids),clipped_sa_read_ids=",".join(sorted(sa_ids)),top_endpoint_peak=best["peak_coordinate"] if best else "NA",
            top_peak_offset=best["peak_coordinate"]-bp if best else "NA",top_peak_unique_reads=best["unique_reads"] if best else 0,
            nearest_n_gap=f"{nearest[0]}-{nearest[1]}" if nearest else "NONE",nearest_n_gap_length=nearest[2] if nearest else 0,
            distance_to_nearest_n_gap=distance(bp,nearest[0],nearest[1]) if nearest else "NA",reads_spanning_any_n_gap=len(gap_union)))
    for rid in sorted(clipped_names):
        for r in sorted(allrecs[rid],key=lambda x:(x.reference_name or "",x.reference_start,x.flag)):
            if r.is_secondary and not a.include_secondary: continue
            qa,qb=qinterval(r); lc,rc=clips(r)
            alignments.append(dict(reference=ds.reference,read_id=rid,alignment_type=atype(r),chrom=r.reference_name,ref_start=r.reference_start+1,
                ref_end=r.reference_end,strand=strand(r),query_start=qa,query_end=qb,mapq=r.mapping_quality,cigar=r.cigarstring,
                left_clip_bp=lc,right_clip_bp=rc,sa_tag=r.get_tag("SA") if r.has_tag("SA") else ""))
    pg="; ".join(" ".join(str(x.get(k,"")) for k in ("ID","PN","VN","CL")).strip() for x in bam.header.to_dict().get("PG",[])) or "Not recorded"
    bam.close(); fa.close()
    return dict(dataset=ds,program=pg,summary=summary,clipped=clipped,sa=sa_rows,endpoint=endpoint,gaps=gaps,alignments=alignments)

def write_tsv(path,rows,fields):
    with open(path,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",extrasaction="ignore"); w.writeheader(); w.writerows(rows)

def main():
    p=argparse.ArgumentParser(description="Independently validate both qHTSF4.1 breakpoints")
    p.add_argument("--config",required=True); p.add_argument("--output-prefix",required=True)
    p.add_argument("--window",type=int,default=10000); p.add_argument("--min-mapq",type=int,default=20)
    p.add_argument("--min-anchor",type=int,default=500); p.add_argument("--span-anchor",type=int,default=500)
    p.add_argument("--min-clip",type=int,default=200); p.add_argument("--cluster-tolerance",type=int,default=50)
    p.add_argument("--min-n-run",type=int,default=50); p.add_argument("--include-secondary",action="store_true")
    a=p.parse_args()
    try: results=[analyse(ds,a) for ds in load_config(a.config)]
    except Exception as e: sys.exit(f"ERROR: {e}")
    prefix=Path(a.output_prefix); prefix.parent.mkdir(parents=True,exist_ok=True)
    keys=["summary","clipped","sa","endpoint","gaps","alignments"]
    groups={k:[row for result in results for row in result[k]] for k in keys}
    fields={
      "summary":["reference","read_sample","chrom","breakpoint_label","configured_breakpoint","continuous_spanning_reads","continuous_read_ids","clipped_reads","clipped_read_ids","clipped_reads_with_sa","clipped_sa_read_ids","top_endpoint_peak","top_peak_offset","top_peak_unique_reads","nearest_n_gap","nearest_n_gap_length","distance_to_nearest_n_gap","reads_spanning_any_n_gap"],
      "clipped":["reference","breakpoint_label","configured_breakpoint","read_id","alignment_type","alignment_location","alignment_strand","endpoint_side","endpoint_coordinate","offset_from_configured_bp","clip_bp","mapq","cigar","query_interval","sa_entry_count"],
      "sa":["reference","breakpoint_label","configured_breakpoint","read_id","source_location","source_type","source_strand","source_cigar","sa_chrom","sa_start","sa_end","sa_strand","sa_cigar","sa_mapq","sa_nm","sa_near_boundary"],
      "endpoint":["reference","breakpoint_label","configured_breakpoint","cluster_rank","cluster_start","cluster_end","peak_coordinate","peak_offset","unique_reads","event_count","read_ids"],
      "gaps":["reference","breakpoint_label","configured_breakpoint","gap_start","gap_end","gap_length","distance_from_breakpoint","breakpoint_inside_gap","gap_spanning_read_count","gap_spanning_read_ids"],
      "alignments":["reference","read_id","alignment_type","chrom","ref_start","ref_end","strand","query_start","query_end","mapq","cigar","left_clip_bp","right_clip_bp","sa_tag"]}
    suffix=dict(summary="breakpoint_summary",clipped="clipped_reads",sa="sa_partners",endpoint="endpoint_clusters",gaps="n_gaps",alignments="read_alignments")
    written=[]
    for k in keys:
        path=f"{prefix}.{suffix[k]}.tsv"; write_tsv(path,groups[k],fields[k]); written.append(path)
    report_path=f"{prefix}.report.md"
    with open(report_path,"w") as o:
        o.write("# qHTSF4.1 independent long-read breakpoint report\n\nEach configured boundary was analysed independently. Counts are unique read counts.\n\n")
        o.write("| Reference | Boundary | Coordinate | Continuous | Clipped | Clipped + SA | Endpoint peak | Peak support | Nearest N gap | Gap distance | N-gap spanning |\n|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|\n")
        for r in groups["summary"]:
            o.write(f"| {r['reference']} | {r['breakpoint_label']} | {r['configured_breakpoint']} | {r['continuous_spanning_reads']} | {r['clipped_reads']} | {r['clipped_reads_with_sa']} | {r['top_endpoint_peak']} | {r['top_peak_unique_reads']} | {r['nearest_n_gap']} | {r['distance_to_nearest_n_gap']} | {r['reads_spanning_any_n_gap']} |\n")
        o.write("\n## Interpretation rules\n\n- Continuous: one alignment spans the configured coordinate with the required anchor on both sides.\n- Clipped + SA: a clipped alignment near the boundary has at least one SA-tag partner.\n- Endpoint peak: the highest-support cluster of clipped alignment starts/ends; it is an estimate, not an automatically confirmed breakpoint.\n- N-gap spanning: an alignment covers both flanks of an N-run; inspect CIGAR and raw-read sequence before treating it as gap-independent evidence.\n\n## BAM alignment provenance\n\n")
        for x in results: o.write(f"- **{x['dataset'].reference}:** `{x['program']}`\n")
    for path in [report_path]+written: print(path)

if __name__=="__main__": main()

