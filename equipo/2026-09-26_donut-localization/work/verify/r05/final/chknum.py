import re,json,math
d=json.load(open('data/paper_numbers.json'))
t=open('paper/generated/numbers.tex',encoding='utf-8').read()
m=dict(re.findall(r'csname (pnum(?:se)?@[^ ]+?).endcsname[{]([^}]*)[}]',t))
def parse(s):
    s=s.replace('$','').replace('{','').replace('}','').replace('\times','x').replace('\text','').replace('−','-').replace('\ensuremath','')
    mm=re.match(r'\s*(-?[\d.]+)\s*(?:x\s*10\^(-?\d+))?',s)
    txt=mm.group(1); e=int(mm.group(2) or 0)
    dec=len(txt.split('.')[1]) if '.' in txt else 0
    return float(txt)*10**e, 10**(e-dec)
bad=0;n=0
for k,v in m.items():
    kind,key=k.split('@',1)
    if key not in d: print('MISSING',key);bad+=1;continue
    ref=d[key]['value'] if kind=='pnum' else d[key].get('se')
    if not isinstance(ref,(int,float)): continue
    val,ulp=parse(v); n+=1
    if abs(val-ref)>0.5*ulp*1.0001: print('MISMATCH',k,v,ref);bad+=1
    if kind=='pnumse':
        pv,pu=parse(m['pnum@'+key])
        sig=len(re.sub(r'^[0.]+','',re.sub(r'x.*','',v.replace('.','')).lstrip('0'))) 
        if abs(pu-ulp)>1e-15*max(pu,ulp) and abs(pu/ulp-1)>1e-6: print('DECIMAL-MISMATCH',key,m['pnum@'+key],v);bad+=1
        digits=round(val/ulp)
        if digits>=25 or digits<1: print('SE-DIGITS',key,v,digits);bad+=1
print('checked',n,'bad',bad)
