#!/shared/apps/sage-7.4/local/bin/sage -python

from sage.all_cmdline import *;

import sys,json;#,linecache,traceback,mongojoin;#os,tempfile,time,datetime,subprocess,signal
from mongojoin.parse import pythonlist2mathematicalist as py2mat;
from mongojoin.parse import mathematicalist2pythonlist as mat2py;
import mongojoin.tools as tools;

#Misc. function definitions
#def PrintException():
#    "If an exception is raised, print traceback of it to output log."
#    exc_type, exc_obj, tb = sys.exc_info();
#    f = tb.tb_frame;
#    lineno = tb.tb_lineno;
#    filename = f.f_code.co_filename;
#    linecache.checkcache(filename);
#    line = linecache.getline(filename, lineno, f.f_globals);
#    print 'EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj);
#    print "More info: ",traceback.format_exc();

#def default_sigpipe():
#    signal.signal(signal.SIGPIPE,signal.SIG_DFL);

'''
def timestamp2unit(timestamp,unit="seconds"):
    if timestamp=="infinite":
        return timestamp;
    else:
        days=0;
        if "-" in timestamp:
            daysstr,timestamp=timestamp.split("-");
            days=int(daysstr);
        hours,minutes,seconds=[int(x) for x in timestamp.split(":")];
        hours+=days*24;
        minutes+=hours*60;
        seconds+=minutes*60;
        if unit=="seconds":
            return seconds;
        elif unit=="minutes":
            return float(seconds)/60.;
        elif unit=="hours":
            return float(seconds)/(60.*60.);
        elif unit=="days":
            return float(seconds)/(60.*60.*24.);
        else:
            return 0;
'''

#def timeleft(starttime,timelimit):
#    "Determine if runtime limit has been reached."
#    if timelimit=="infinite":
#        return True;
#    else:
#        return ((timestamp2unit(timelimit)-(time.time()-starttime)))>0;

#docsfile=sys.argv[1];
#basecoll=sys.argv[2];
#dbindexes=sys.argv[3:];
#with open(docsfile,'r') as docstream:
#for line in sys.stdin:#docstream:
#for line in iter(sys.stdin.readline, ''):
#line=sys.stdin.readline();
dbtype=sys.argv[1];
dbusername=sys.argv[2];
dbpassword=sys.argv[3];
dbhost=sys.argv[4];
dbport=sys.argv[5];
dbname=sys.argv[6];
if dbtype=="mongodb":
    import mongojoin;
    if dbusername==None:
        dbclient=mongojoin.MongoClient("mongodb://"+dbhost+":"+dbport+"/"+dbname);
    else:
        dbclient=mongojoin.MongoClient("mongodb://"+dbusername+":"+dbpassword+"@"+dbhost+":"+dbport+"/"+dbname+"?authMechanism=SCRAM-SHA-1");
    db=dbclient[dbname];

for line in iter(sys.stdin.readline,''):
    #with open("/gss_gpfs_scratch/altman.ro/SLURMongo/scripts/1.err","w") as errstream:
    #    errstream.write(line);
    #    errstream.flush();
    polydoc=json.loads(line.rstrip("\n"));

    #with open(workpath+"/"+jobstepname+".docs","r") as docstream:
    #    count=0;
    #    #polydoc=next(polycurs,None);
    #    line=docstream.readline();
    #    while line!="":# and timeleft(starttime,timelimit):
    #        polydoc=json.loads(line.rstrip("\n"));
    #        #iddoc=dict([(x,polydoc[x]) for x in dbindexes if x in polydoc.keys()]);
    #Read in pertinent fields from JSON
    polyid=polydoc['POLYID'];
    nverts=mat2py(polydoc['NVERTS']);

    lp=LatticePolytope(nverts);
    dlp=LatticePolytope(lp.polar().normal_form());

    dverts=[list(x) for x in dlp.vertices().column_matrix().columns()];

    #lp_facets=lp.faces_lp(codim=1);

    #lp_facetbndrypts=[[list(y) for y in x.boundary_points()] for x in lp_facets];

    #lp_bndrypts_dup=[y for x in lp_facetbndrypts for y in x];
    #lp_bndrypts=[lp_bndrypts_dup[i] for i in range(len(lp_bndrypts_dup)) if lp_bndrypts_dup[i] not in lp_bndrypts_dup[:i]];

    #dlp_facets=dlp.faces_lp(dim=3);

    #dlp_facetbndrypts=[[list(y) for y in x.boundary_points()] for x in dlp_facets];
    #dlp_facetinterpts=[[list(y) for y in x.interior_points()] for x in dlp_facets];
    #dlp_maxcone_normalform_dup=[[list(y) for y in LatticePolytope(x.vertices().column_matrix().columns()+[vector((0,0,0,0))]).normal_form().column_matrix().columns()] for x in dlp_facets];
    #dlp_maxcone_normalform_inds=[i for i in range(len(dlp_maxcone_normalform_dup)) if dlp_maxcone_normalform_dup[i] not in dlp_maxcone_normalform_dup[:i]];
    
    #dlp_bndrypts_dup=[y for x in dlp_facetbndrypts for y in x];
    #dlp_bndrypts=[dlp_bndrypts_dup[i] for i in range(len(dlp_bndrypts_dup)) if dlp_bndrypts_dup[i] not in dlp_bndrypts_dup[:i]];

    #dlp_interpts_dup=[y for x in dlp_facetinterpts for y in x];
    #dlp_interpts=[dlp_interpts_dup[i] for i in range(len(dlp_interpts_dup)) if dlp_interpts_dup[i] not in dlp_interpts_dup[:i]];

    nformlist=[];
    nform2skellist=[]
    nformcountlist=[];
    faceinfolist=[];
    for facet in dlp.faces_lp(dim=3):
        nform=[list(w) for w in LatticePolytope(facet.vertices().column_matrix().columns()+[vector((0,0,0,0))]).normal_form().column_matrix().columns() if w!=vector((0,0,0,0))];
        if nform not in nformlist:
            nform2skel=[list(w) for w in LatticePolytope(nform).boundary_points().column_matrix().columns()];
            nformlist+=[nform];
            nform2skellist+=[nform2skel];
            nformcountlist+=[1];
            dim0=[facet.nvertices()];
            face1intpts=[len(y.interior_points()) for y in facet.faces_lp(dim=1)];
            face1bdrypts=[len(y.boundary_points()) for y in facet.faces_lp(dim=1)];
            dim1int=[len(face1intpts),sum(face1intpts),min(face1intpts),max(face1intpts)];
            dim1bdry=[len(face1bdrypts),sum(face1bdrypts),min(face1bdrypts),max(face1bdrypts)];
            face2intpts=[len(y.interior_points()) for y in facet.faces_lp(dim=2)];
            face2bdrypts=[len(y.boundary_points()) for y in facet.faces_lp(dim=2)];
            dim2int=[len(face2intpts),sum(face2intpts),min(face2intpts),max(face2intpts)];
            dim2bdry=[len(face2bdrypts),sum(face2bdrypts),min(face2bdrypts),max(face2bdrypts)];
            faceinfolist+=[dim0+dim1int+dim1bdry+dim2int+dim2bdry];
        else:
            nformcountlist[nformlist.index(nform)]+=1;
    
    nformexistlist=[x["NFORM"] for x in list(db["FACET"].find({"NFORM":{"$in":[py2mat(nform) for nform in nformlist]}},{"_id":0,"NFORM":1}).limit(len(nformlist)))];
    #print("+POLY."+json.dumps({'POLYID':polyid},separators=(',',':'))+">"+json.dumps({'DVERTS':py2mat(dverts)},separators=(',',':')));
    #sys.stdout.flush();
    
    facetlist=[];
    for i in range(len(nformlist)):
        nform=nformlist[i];
        nform2skel=nform2skellist[i];
        faceinfo=faceinfolist[i];
        ninstances=nformcountlist[i];
        facetlist+=[{'NFORM':py2mat(nform),'NINST':ninstances}];
        #print("&FACET."+json.dumps({'NFORM':py2mat(nform)},separators=(',',':'))+">"+json.dumps({'POS':{'POLYID':polyid,'NINST':ninstances}},separators=(',',':')));
        mat_nform=py2mat(nform);
        if mat_nform not in nformexistlist:
            print("+FACET."+json.dumps({'NFORM':mat_nform},separators=(',',':'))+">"+json.dumps({'NFORM2SKEL':py2mat(nform2skel),'FACEINFO':py2mat(faceinfo),'facetfinetriangsMARK':False,'facetallfinetriangsMARK':False},separators=(',',':')));
            sys.stdout.flush();

    print("+POLY."+json.dumps({'POLYID':polyid},separators=(',',':'))+">"+json.dumps({'DVERTS':py2mat(dverts),'FACETLIST':facetlist},separators=(',',':')));
    print("@");#+basecoll+"."+json.dumps(dict([(x,polydoc[x]) for x in dbindexes]),separators=(',',':')));
    sys.stdout.flush();
    #line=sys.stdin.readline();

#        line=docstream.readline();
dbclient.close();
#if __name__ == "__main__":
#    main();

#with tempfile.NamedTemporaryFile(dir=workpath,delete=False) as tempstream:
#    while line!="":
#        tempstream.write(line);
#        tempstream.flush();
#        line=docstream.readline();
#    os.rename(tempstream.name,docstream.name);

#if line!="" and not timeleft(starttime,timelimit):   
#submit=subprocess.Popen("sbatch "+mainpath+"/"+jobnumpad+"/maxcones_"+jobnumpad+".job",shell=True,stdout=subprocess.PIPE,preexec_fn=default_sigpipe);
#submitcomm=submit.communicate()[0].rstrip("\n");
#Print information about job submission
#print "";
#print datetime.datetime.now().strftime("%Y %m %d %H:%M:%S");
#print "Res"+submitcomm[1:]+" as maxcones_"+jobstepname+".";
#print "";
#print "";
#print "Paused"
#sys.stdout.flush();

#docstream.close();
#except Exception as e:
#    PrintException();