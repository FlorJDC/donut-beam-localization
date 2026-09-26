import fitz, glob, os
out='png'
for f in sorted(glob.glob('clone/paper/figures/fig*.pdf')):
    d=fitz.open(f); p=d[0]
    print(os.path.basename(f), p.rect.width/72*2.54, 'cm wide', p.rect.height/72*2.54,'cm high')
    p.get_pixmap(dpi=130).save(os.path.join(out, os.path.basename(f)[:-4]+'.png'))
d=fitz.open('clone/paper/main.pdf'); print('main pages', len(d))
for i,p in enumerate(d):
    imgs=p.get_images(); dr=len(p.get_drawings()) if i<0 else 0
    p.get_pixmap(dpi=110).save(os.path.join(out,'main_p%02d.png'%(i+1)))
