import sys, json, math, re, ezdxf
from shapely.geometry import Polygon, MultiPolygon, Point, mapping, shape
from shapely.ops import unary_union
from shapely import affinity
dxf, immjs, out = sys.argv[1], sys.argv[2], sys.argv[3]
doc = ezdxf.readfile(dxf); msp = doc.modelspace()

# ---- muri: campiture dei blocchi Muro/Forometria
polys = []
openings = []
for e in msp.query('INSERT[layer=="Ferrandes_Muri"]'):
    for v in e.virtual_entities():
        if v.dxftype()=='HATCH':
            for p in v.paths:
                vs = getattr(p,'vertices',None)
                if vs and len(vs)>=3:
                    pg = Polygon([(x[0],x[1]) for x in vs]).buffer(0)
                    if pg.area>1e-6: polys.append(pg)
        elif v.dxftype()=='INSERT' and v.dxf.layer in ('Porte Archicad','Finestre Archicad'):
            xs=[];ys=[]
            for y in v.virtual_entities():
                if y.dxftype()=='LINE':
                    for q in (y.dxf.start,y.dxf.end): xs.append(q[0]); ys.append(q[1])
            ins=v.dxf.insert; rot=math.radians(v.dxf.rotation)
            # larghezza: estensione lungo l'asse del blocco
            ux,uy=math.cos(rot),math.sin(rot)
            proj=[(x-ins[0])*ux+(y-ins[1])*uy for x,y in zip(xs,ys)]
            w = (max(proj)-min(proj)) if proj else 1
            mid = (min(proj)+max(proj))/2 if proj else 0
            openings.append({'kind':'porta' if v.dxf.layer=='Porte Archicad' else 'finestra','name':v.dxf.name,
                             'c':(ins[0]+ux*mid, ins[1]+uy*mid),'w':w,'u':(ux,uy)})
walls = unary_union(polys)
print('muri', walls.geom_type, round(walls.area,1), 'm2; aperture', len(openings))

# ---- vani
rooms=[]
labels=[(e.dxf.insert, dict((a.dxf.tag,a.dxf.text) for a in e.attribs).get('ROOM_NAME',''),
         [re.sub(r'\\[A-Za-z][^;]*;|[{}]|\\S[^;]*;','',y.text).strip() for y in doc.blocks[e.dxf.name] if y.dxftype()=='MTEXT'])
        for e in msp.query('INSERT[layer=="Ferrandes_Zona"]')]
for e in msp.query('LWPOLYLINE[layer=="Ferrandes_Zona"]'):
    pg=Polygon([p[:2] for p in e.get_points()]).buffer(0)
    lab=[l for l in labels if pg.contains(Point(l[0][0],l[0][1]))]
    nm, area = (lab[0][1], lab[0][2][0].split(' ')[0]) if lab else ('', '')
    rooms.append({'poly':pg,'nome':nm,'mq':area.replace(',','.')})

# ---- varchi: muro "chiuso" meno muro, vicino a ogni apertura
gaps=[]
for o in openings:
    r = o['w']/2+.05
    zone = Point(o['c']).buffer(r+.4)
    closed = walls.intersection(zone).buffer(r, join_style=2).buffer(-r, join_style=2)
    g = closed.difference(walls).intersection(Point(o['c']).buffer(r+.15))
    g = g.buffer(.02).buffer(-.02)
    if not g.is_empty and g.area>.05: gaps.append((o,g))
print('varchi trovati', len(gaps), '/', len(openings))

# ---- sagoma catastale (metri locali)
s=open(immjs).read(); im=[i for i in json.loads(s[s.index('['):].rstrip().rstrip(';')) if i['id']=='66_224'][0]
lng0,lat0=im['centro']; kx=111320*math.cos(math.radians(lat0)); ky=110540
cad = unary_union([Polygon([((c[0]-lng0)*kx,(c[1]-lat0)*ky) for c in poly[0]]) for poly in im['sagoma']['coordinates']]).buffer(0)
# involucro = muri + vani + varchi (i soli muri non racchiudono i vani: le porte li lasciano aperti)
env = unary_union([walls]+[r['poly'] for r in rooms]+[g for _,g in gaps]).buffer(.25, join_style=2).buffer(-.25, join_style=2)
env = unary_union([Polygon(g.exterior) for g in (env.geoms if hasattr(env,'geoms') else [env])])
print('sup. DXF', round(env.area,1), 'catasto', round(cad.area,1))

def place(g, th, tx, ty): return affinity.translate(affinity.rotate(g, th, origin=(0,0)), tx, ty)
c0 = env.centroid; base = affinity.translate(env, -c0.x, -c0.y)
cc = cad.centroid
def iou(th,tx,ty):
    g=place(base,th,tx,ty); return g.intersection(cad).area/g.union(cad).area
best=max(((iou(th,cc.x,cc.y),th) for th in range(0,360,2)))
th=best[1]; tx,ty=cc.x,cc.y; sc=best[0]
for step in (1,.3,.1,.03):
    improved=True
    while improved:
        improved=False
        for d in [(step,0,0),(-step,0,0),(0,step,0),(0,-step,0),(0,0,step),(0,0,-step)]:
            v=iou(th+d[0]*3,tx+d[1],ty+d[2])
            if v>sc+1e-6: sc,th,tx,ty=v,th+d[0]*3,tx+d[1],ty+d[2]; improved=True
# prova anche speculare
print('allineamento: rotazione %.1f° IoU %.3f'%(th,sc))
T=lambda g: place(affinity.translate(g,-c0.x,-c0.y),th,tx,ty)
G=lambda g: shape(json.loads(json.dumps(mapping(affinity.affine_transform(g,[1/kx,0,0,1/ky,lng0,lat0])))))

def feat(g, props):
    g=G(T(g)); return {'type':'Feature','properties':props,'geometry':json.loads(json.dumps(mapping(g)))}
F=[feat(walls,{'k':'muri'}), feat(env,{'k':'pavimento'})]
for o,g in gaps: F.append(feat(g,{'k':o['kind'],'w':round(o['w'],2)}))
# pavimento netto: tutto l'interno dell'involucro tranne i muri (vani, soglie di porte e finestre, passaggi)
net = env.difference(walls).buffer(.03, join_style=2).simplify(.005)
for g in (net.geoms if hasattr(net,'geoms') else [net]):
    if g.area > .02: F.append(feat(g,{'k':'pavimento_netto'}))
print('pavimento netto', round(net.area,1), 'm2 in', len(getattr(net,'geoms',[net])), 'parti')
for r in rooms:
    lp=r['poly'].representative_point()
    F.append(feat(r['poly'],{'k':'vano','nome':r['nome'],'mq':r['mq']}))
res={'type':'FeatureCollection','features':F,'meta':{'fonte':'Ferrandes F66-P224 stato di fatto.dxf','rotazione':round(th,2),'iou':round(sc,3)}}
open(out,'w').write('window.PLANIMETRIE = window.PLANIMETRIE || {};\nwindow.PLANIMETRIE["66_224"] = '+json.dumps(res,separators=(',',':'))+';\n')
print('scritto', out)
