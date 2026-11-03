from pathlib import Path
import urllib.request,json,hashlib,re
ROOT=Path(__file__).resolve().parents[1]
REPO='jbloomlab/Ab-CGGnaive_DMS'; COMMIT='5ffc88950f760a297e2196ac449111e3db023e80'
PATHS=['README.md','config.yaml','Snakefile','compute_binding_Kd.Rmd','compute_expression_meanF.Rmd','collapse_scores.Rmd','data/wildtype_sequence.fasta','data/CGGnaive_sites.csv','data/barcode_runs.csv','results/summary/compute_expression_meanF.md','results/summary/Titeseq-modeling.md','results/binding_Kd/bc_binding.csv','results/binding_Kd/variant_binding.csv','results/expression_meanF/bc_expression.csv','results/expression_meanF/variant_expression.csv','results/variants/codon_variant_table.csv','results/final_variant_scores/final_variant_scores.csv']
def get(url):
 with urllib.request.urlopen(url,timeout=90) as r:return r.read()
def blobsha(b):return hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
def main():
 tree=json.loads(get(f'https://api.github.com/repos/{REPO}/git/trees/{COMMIT}?recursive=1'))
 assert not tree.get('truncated')
 blobs={e['path']:e for e in tree['tree'] if e['type']=='blob'}
 records=[]
 for name in PATHS:
  raw=get(f'https://raw.githubusercontent.com/{REPO}/{COMMIT}/{name}')
  if blobsha(raw)!=blobs[name]['sha']:raise RuntimeError('Git mismatch '+name)
  pointer=raw.startswith(b'version https://git-lfs.github.com/spec/v1')
  b=raw
  if pointer:
   oid=re.search(rb'oid sha256:([0-9a-f]+)',raw).group(1).decode(); size=int(re.search(rb'size (\d+)',raw).group(1))
   b=get(f'https://media.githubusercontent.com/media/{REPO}/{COMMIT}/{name}')
   if len(b)!=size or hashlib.sha256(b).hexdigest()!=oid:raise RuntimeError('LFS mismatch '+name)
  p=ROOT/'inputs'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
  records.append({'repository':REPO,'commit':COMMIT,'path':name,'git_blob':blobs[name]['sha'],'lfs':pointer,'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b)})
  print(name,len(b),flush=True)
 (ROOT/'results/SOURCE_MANIFEST.json').write_text(json.dumps({'files':records},indent=2)+'\n')
if __name__=='__main__':main()
