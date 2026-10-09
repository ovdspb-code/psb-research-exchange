"""New data extraction and identifiability/precision calculations; never executes R1."""
from pathlib import Path
from itertools import product
from collections import Counter,defaultdict
from fractions import Fraction
from math import ceil,log,sqrt,exp,isfinite
import csv,json,hashlib,statistics,re,sys
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'source'; IN=ROOT/'inputs'; OUT=ROOT/'results'
CODONS=[''.join(x) for x in product('TCAG',repeat=3)]
AAS='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
GENETIC=dict(zip(CODONS,AAS))
def translate(nt):
 nt=nt.upper()
 if len(nt)%3:raise ValueError('Frame')
 return ''.join(GENETIC[nt[i:i+3]] for i in range(0,len(nt),3))
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def write(name,rows,fields=None):
 if fields is None:fields=list(rows[0])
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def required(ok,msg):
 if not ok:raise RuntimeError(msg)
def jsonout(name,data):(OUT/name).write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
def finite(v):
 try:return isfinite(float(v))
 except (ValueError,TypeError):return False

def evidence():
 manifest=json.loads((OUT/'SOURCE_MANIFEST.json').read_text())
 for e in manifest['files']:
  b=(IN/e['path']).read_bytes();required(len(b)==e['bytes'] and hashlib.sha256(b).hexdigest()==e['sha256'],'input bytes')
 starts=read(SOURCE/'planned_starting_constructs.csv')
 focal=[x for x in starts if x['role'] in ['FOCAL_ancestor','FOCAL_current','FOCAL_rewrite']]
 required(len(focal)==3,'focal roles')
 sites=read(IN/'data/CGGnaive_sites.csv')
 # The site table's amino acids are independently checked against the complete assay DNA.
 fasta=''.join((IN/'data/wildtype_sequence.fasta').read_text().splitlines()[1:]).upper(); aa=translate(fasta)
 required(len(aa)==235,'assay boundaries')
 naive=''.join(x['amino_acid'] for x in sorted(sites,key=lambda x:int(x['site_scFv'])))
 required(aa==naive,'site table does not match assay')
 mapping=[]
 for s in focal:
  h=translate(s['hc']);l=translate(s['lc'])
  required(h==aa[:112] and len(l)==107,'H chain / L length')
  diffs=[(i+1,x,y) for i,(x,y) in enumerate(zip(aa[127:234],l)) if x!=y]
  required(diffs==[(36,'Y','N')],'protein mutation')
  mapping.append({'role':s['role'],'length_H':len(h),'length_L':len(l),'paired_protein_id':s['protein_id'],'scFv_position':163,'IMGT':'Y42(L)N','LC_S93_codon':s['lc'][276:279],'variable_AA_same':True})
 write('FOCAL_MAPPING.csv',mapping)
 bcmap=read(IN/'results/variants/codon_variant_table.csv')
 selected_map=[x for x in bcmap if x['aa_substitutions']=='Y163N']
 lookup={(x['library'],x['barcode']):x for x in selected_map}
 bcsummary=[]
 for filename,value in [('results/binding_Kd/bc_binding.csv','log10Ka'),('results/expression_meanF/bc_expression.csv','expression')]:
  rows=[x for x in read(IN/filename) if x['aa_substitutions']=='Y163N' or x['variant_class']=='wildtype']
  write(Path(filename).stem+'_Y163N_WT.csv',rows)
  for lib in ['lib1','lib2']:
   for kind in ['Y163N','WT']:
    sub=[x for x in rows if x['library']==lib and ((x['aa_substitutions']=='Y163N') if kind=='Y163N' else x['variant_class']=='wildtype')]
    vals=[float(x[value]) for x in sub if finite(x[value])]
    bcsummary.append({'file':filename,'library':lib,'variant':kind,'rows':len(sub),'finite_fits':len(vals),'missing_fits':len(sub)-len(vals),'mean_legacy_fit':statistics.mean(vals) if vals else None,'interpretation':'legacy individual fits; not replicate of current pooled-variant fit; no biological CI'})
 write('BARCODE_FIT_QA.csv',bcsummary)
 write('Y163N_BARCODE_CODONS.csv',selected_map)
 vbind=[x for x in read(IN/'results/binding_Kd/variant_binding.csv') if x['aa_substitutions'] in ['Y163N','WT']]
 vexpr=[x for x in read(IN/'results/expression_meanF/variant_expression.csv') if x['aa_substitutions'] in ['Y163N','WT']]
 library=[]
 for lib in ['lib1','lib2']:
  b={x['aa_substitutions']:x for x in vbind if x['library']==lib};e={x['aa_substitutions']:x for x in vexpr if x['library']==lib}
  library.append({'library':lib,'variant':'Y163N','WT_binding_log10Ka':float(b['WT']['bind_CGG']),'Y163N_binding_log10Ka':float(b['Y163N']['bind_CGG']),'delta_binding_log10Ka':float(b['Y163N']['bind_CGG'])-float(b['WT']['bind_CGG']),'WT_expression_ln':float(e['WT']['expression']),'Y163N_expression_ln':float(e['Y163N']['expression']),'delta_expression_ln':float(e['Y163N']['expression'])-float(e['WT']['expression']),'binding_barcode_support':int(float(b['Y163N']['n_bc_bind_CGG'])),'expression_barcode_support':int(e['Y163N']['n_bc_expr'])})
 write('Y163N_WT_LIBRARY_EVIDENCE.csv',library)
 delta_b=statistics.mean(x['delta_binding_log10Ka'] for x in library);delta_e=statistics.mean(x['delta_expression_ln'] for x in library)
 required(abs(delta_b-.25938)<.000005 and abs(delta_e+1.26992)<.000005,'final score match')
 match_rows=[]
 for s in selected_map:
  seq=list(fasta)
  for old,pos,new in re.findall(r'([ACGT]{3})(\d+)([ACGT]{3})',s['codon_substitutions']):
   i=3*(int(pos)-1);required(''.join(seq[i:i+3])==old,'codon reference');seq[i:i+3]=list(new)
  seq=''.join(seq)
  for f in focal:
   match_rows.append({'library':s['library'],'barcode':s['barcode'],'role':f['role'],'exact_variable_DNA_match':seq[:336]==f['hc'] and seq[381:702]==f['lc']})
 write('EXACT_SPELLING_MATCH.csv',match_rows)
 required(not any(x['exact_variable_DNA_match'] for x in match_rows),'unexpected exact spelling match')
 evid=[{'object':'Y163N / Y42(L)N','evidence':'DIRECT','assay':'yeast scFv: 235 residues','readout':'library-pooled binding and expression','support':'lib1 20 + lib2 24 barcodes','limitation':'two libraries; not 44 biological replicates; no CI asserted'},
 {'object':'a/x/R exact DNA','evidence':'ABSENT','assay':'native receptor / matched exact coding','readout':'binding, surface expression, viability and declared competence','support':'0 exact variable DNA matches among 44 mapped Y163N barcodes','limitation':'same AA does not prove same codon-dependent expression'},
 {'object':'Y36N + future substitution','evidence':'TRANSPORTED','assay':'single-mutant DMS additive extrapolation','readout':'joint usefulness statuses','support':'accepted 002 conditional map','limitation':'mature double-mutant binding/expression not measured by single-mutant row'},
 {'object':'LC mutation supply','evidence':'TRANSPORTED','assay':'mouse nonfunctional kappa repertoire model to focal active receptor','readout':'two-event first protein exit','support':'MK_RS5NF model used in 001/002','limitation':'no actual ordered focal event observations'}]
 write('COMPETENCE_EVIDENCE.csv',evid)
 jsonout('COMPETENCE_RESULTS.json',{'focal':mapping,'assay':{'HC_residues':[1,112],'linker_residues':[113,127],'LC_shared_residues':[128,234],'extra_terminal_residue':aa[234],'linker_AA':aa[112:127]},'libraries':library,'mean_delta_binding_log10Ka':delta_b,'mean_delta_expression_ln':delta_e,'relative_expression_ratio_interpretation':exp(delta_e),'Y163N_barcode_map_rows':len(selected_map),'exact_DNA_matches':sum(x['exact_variable_DNA_match'] for x in match_rows),'raw_FCS_or_sequencing_reanalysed':False,'individual_binding_file_does_not_reproduce_current_pooled_fit':True})

def observation():
 weights=[x for x in read(SOURCE/'all_destination_weights.csv') if x['case_id']=='D20_M23_GC102__seq7__seq6']
 panel={(x['position_chain'],x['target_aa']):int(x['expected_joint_status']) for x in read(SOURCE/'planned_phenotype_panel.csv') if x['role']=='FOCAL_RESERVE'}
 required(len(weights)==2140 and len(panel)==85,'retained scope')
 masses={}
 for role in ['ancestor','current','rewrite']:
  masses[role]={'known_useful':0.,'unknown':0.,'known_nonuseful':0.}
  for x in weights:
   v=float(x['w2_'+role]);status=0 if x['stop']=='True' else panel.get((x['position_chain'],x['target_aa']))
   masses[role]['unknown' if status is None else 'known_useful' if status==1 else 'known_nonuseful']+=v
  masses[role]['no_exit']=1-sum(masses[role].values())
  required(masses[role]['no_exit']>=0,'no exit mass')
 coarse={name:masses[a]['known_useful']-masses[b]['known_useful']-masses[b]['unknown'] for name,a,b in [('loss','ancestor','current'),('rewrite','rewrite','current'),('restore','rewrite','ancestor')]}
 required(all(x<0 for x in coarse.values()),'unexpected coarse feasibility')
 write('COARSE_READOUT_MODEL_MASS.csv',[{'role':r,**v} for r,v in masses.items()])
 # Same source codon at focal S93, all exactly-two single-base paths in this codon.
 # Includes reversals; does not simulate SHM or assert biological frequencies.
 paths=[]
 for role,start in [('ancestor','AGC'),('current','AGT'),('rewrite','AGC')]:
  for first in CODONS:
   if sum(a!=b for a,b in zip(start,first))!=1:continue
   for final in CODONS:
    if sum(a!=b for a,b in zip(first,final))!=1:continue
    first_exit=GENETIC[first] if GENETIC[first]!=GENETIC[start] else GENETIC[final] if GENETIC[final]!=GENETIC[start] else 'NO_EXIT'
    paths.append({'role':role,'start':start,'intermediate':first,'endpoint':final,'first_protein_exit':first_exit,'endpoint_AA':GENETIC[final],'endpoint_hamming':sum(a!=b for a,b in zip(start,final))})
 write('S93_EXACT_TWO_EVENT_PATHS.csv',paths)
 groups=defaultdict(set)
 for x in paths:groups[(x['role'],x['endpoint'])].add(x['first_protein_exit'])
 ambiguous=[{'role':r,'endpoint':e,'possible_first_exits':'|'.join(sorted(exits))} for (r,e),exits in groups.items() if len(exits)>1]
 write('ENDPOINT_AMBIGUITIES.csv',ambiguous)
 witnesses=[['AGT','AAT','AAA'],['AGT','AGA','AAA']]
 required([translate(p[1]) for p in witnesses]==['N','R'],'exit control')
 required(all(p[-1]=='AAA' for p in witnesses),'endpoints')
 # Independent rational enumeration of a finite synthetic chain vs matrix recurrence.
 # s0 -> useful .1, failure .2, silent s1 .7; s1 -> useful .3, failure .4, silent s0 .3.
 b=[Fraction(1,10),Fraction(3,10)]; T=[[0,Fraction(7,10)],[Fraction(3,10),0]]
 rec=b[0]+sum(T[0][j]*b[j] for j in range(2));enum=Fraction(1,10)+Fraction(7,10)*Fraction(3,10)
 required(rec==enum==Fraction(31,100),'synthetic readout recurrence')
 jsonout('OBSERVATION_RESULTS.json',{'status':'ENDPOINT_ONLY_NO_GO; ORDERED_OR_VALIDATED_KERNEL_CONDITIONAL','exact_paths':len(paths),'ambiguous_endpoint_groups':len(ambiguous),'focal_S93_source_codons':{'ancestor':'AGC','current':'AGT','rewrite':'AGC'},'witnesses':witnesses,'synthetic_enum_phi2':str(enum),'synthetic_scope':'instrument check only; not biological observations','model_only_group_masses':masses,'model_only_coarse_contrast_lower':coarse,'coarse_85_status_interface_sufficient':False,'coarse_scope':'model illustrates loss of common-label cancellation; not measured probabilities'})

def precision():
 k=2141;eps=.00095;alpha=.025
 rows=[]
 # Known margins conditional on the old model AND the as-yet-unmeasured planned labels.
 margins={'loss':.00192157,'rewrite':.00806796,'restore':.00232927}
 for roles in [3,7]:
  for bias in [0,.0001,.0003,.0005,.00095,.001]:
   stat=eps-bias
   for method in ['subset_Hoeffding_TV','expected_L1_McDiarmid_TV']:
    if stat<=0:n=None
    elif method=='subset_Hoeffding_TV':n=ceil((k*log(2)+log(roles/alpha))/(2*stat**2))
    else:n=ceil((sqrt(k-1)/2+sqrt(log(roles/alpha)/2))**2/stat**2)
    rows.append({'interface':method,'roles':roles,'alpha_S':alpha,'categories':k,'total_TV_radius':eps,'bias_budget_per_role':bias,'stat_radius':max(0,stat),'sufficient_iid_histories_per_role':n,'independent_units':'histories independent under fixed start/carrier law; NOT reads or clonal descendants','meaning':'sufficient bound only; at or above bias radius sampling cannot certify','additional_assumption':'certified TV bias bound covers misclassification+censoring+transfer+batch'})
 # Full useful labels scenario is NOT supplied by 85-panel; avoids dimension only with added phenotypes.
 for name,m in margins.items():
  n=ceil(2*log(6/alpha)/m**2)
  rows.append({'interface':'ALL_LABELS_KNOWN_Bernoulli_Hoeffding_HYPOTHETICAL','roles':3,'alpha_S':alpha,'categories':2,'total_TV_radius':m/2,'bias_budget_per_role':0,'stat_radius':m/2,'sufficient_iid_histories_per_role':n,'independent_units':'independent ordered histories','meaning':'contrast '+name+'; replace m by observed positive margin; strict signs require slightly larger n','additional_assumption':'every actual first-exit useful status independently known; NOT the frozen 85-panel'})
 write('PRECISION_FEASIBILITY.csv',rows)
 by_n=[]
 for n in [10000,100000,1000000,10000000,100000000,1000000000]:
  tv=sqrt(k-1)/(2*sqrt(n))+sqrt(log(3/alpha)/(2*n))
  by_n.append({'n_independent_per_focal_start':n,'sufficient_distribution_free_TV_radius':min(1.,tv),'contrast_uncertainty_from_supply':min(2.,2*tv),'if_all_labels_known_Bernoulli_radius':sqrt(log(6/alpha)/(2*n))})
 write('PRECISION_BY_N.csv',by_n)
 # Synthetic confidence coverage check: Hoeffding union bounds algebra is checked, not simulated.
 required(6*exp(-2*1000000*(sqrt(log(6/alpha)/(2*1000000)))**2)<=alpha+1e-12,'union bound')
 jsonout('PRECISION_RESULTS.json',{'epsilon':eps,'minimum_conditional_margin':min(margins.values()),'residual_smallest_margin_after_2epsilon':min(margins.values())-2*eps,'per_role_sufficient_n_iid_TV':rows[0]['sufficient_iid_histories_per_role'],'sample_count_is_not_lower_bound':True,'event_capture_efficiency_sensitivity':{str(p):ceil(rows[0]['sufficient_iid_histories_per_role']/p) for p in [1,.5,.1,.01]},'efficiency_scope':'conditional mean screened-history planning only under capture independent of outcome; no high-probability collection guarantee; no known physical capture rate','cluster_note':'n counts independent sampled lineages/batches under declared law; repeated reads reduce measurement error only; an n_eff heuristic is not a certificate'})

def main():
 OUT.mkdir(exist_ok=True);evidence();observation();precision()
 jsonout('RUN_SUCCESS.json',{'status':'FEASIBILITY_CALCULATIONS_COMPLETE','exit':0,'old_runner_executed':False,'biological_measurements_generated':False,'python':sys.version,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
 print(json.dumps({'status':'FEASIBILITY_CALCULATIONS_COMPLETE','matched_protein':'Y163N','biological_GO':False,'endpoint_only':'NO_GO','selected_consultation_route':'exact-DNA competence + selection-decoupled ordered-event or validated transition readout','R1_modified':False}))
if __name__=='__main__':
 (OUT/'RUN_SUCCESS.json').unlink(missing_ok=True)
 main()
