#  Structure of the dti package  version 1.1-5

#
#   Class-definitions  in classes.r
#

#
#   R methods and functs 
#
Rstructure <- list(
#
#   File io.r
#
list(funct="dtiData",infile="io.r",visible=TRUE,callsR=c("replind","create.designmatrix.dti"),callsF=c("initdata"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="readDWIdata",infile="io.r",visible=TRUE,callsR=c("replind","create.designmatrix.dti","vcrossp"),callsF=c("initdata"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="print",infile="io.r",visible=TRUE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="show",infile="io.r",visible=TRUE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="summary",infile="io.r",visible=TRUE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="tensor2medinria",infile="io.r",visible=TRUE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="medinria2tensor",infile="io.r",visible=TRUE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
#
#     File dti.r
#
list(funct="dtiTensor",infile="dti.r",visible=TRUE,callsR=c("sdpar","sioutlier","connect.mask","dti3Dreg","plmatrix","pnlrdtirg","opttensR","rho2D","tensRres","pmatrix","pnltens","mcorr","dti3Dev"),callsF=c("nlrdtirg"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="opttensR",infile="dti.r",visible=FALSE,callsR=c(NULL),callsF=c("opttensR"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="tensRres",infile="dti.r",visible=FALSE,callsR=c(NULL),callsF=c("tensRres"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="rho2D",infile="dti.r",visible=FALSE,callsR=c(NULL),callsF=c("rho2D0"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="D2rho",infile="dti.r",visible=FALSE,callsR=c(NULL),callsF=c("D2rho0"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dtiIndices",infile="dti.r",visible=TRUE,callsR=c("dtiind3D"),callsF=c("nlrdtirg"),callsC=c(NULL),purpose="",moveto=""), 
#
#     File hardi.r
#  
list(funct="dwiQball",infile="hardi.r",visible=TRUE,callsR=c("sioutlier","connect.mask","design.spheven","plzero","plzeroaganji","mcorr"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="design.spheven",infile="hardi.r",visible=FALSE,callsR=c("sphcoord","getsphericalharmonicseven"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="plzero",infile="hardi.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="plzeroaganji",infile="hardi.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
#
#     File spherharm.r
#  
list(funct="getsphericalharmonicseven",infile="spherharm.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
#
#    File plot.r
#
list(funct="plot",infile="plot.r",visible=FALSE,callsR=c("identifyFA"),callsF=c("dti2Dga","dti2Dfa"),callsC=c(NULL),purpose="",moveto=""), 
#
#    File show3d.r
#
list(funct="show3d",infile="show3d.r",visible=FALSE,callsR=c("selectCube","extract","show3dData","show3dODF","show3dTens","tracking","expandFibers","design.spheven","dkiIndices","dkiDesign"),callsF=c("mixtradi","mixandir"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="show3dTens",infile="show3d.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="show3dData",infile="show3d.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="show3dODF",infile="show3d.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
#
#    File utilities.r
#
list(funct="sdpar",infile="utilities.r",visible=TRUE,callsR=c("awslinsd"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getsdofsb",infile="utilities.r",visible=FALSE,callsR=c("extract","awslinsd"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="[",infile="utilities.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="combineDWIdata",infile="utilities.r",callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="subsetg",infile="utilities.r",callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="extract",infile="utilities.r",visible=FALSE,callsR=c("dti3Dall","dti3Dev","dti3Dand"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getmask",infile="utilities.r",callsR=c(NULL),callsF=c("getmask"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="selectCube",infile="utilities.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
#
#    File misc.r
#
list(funct="sioutlier",infile="misc.r",visible=FALSE,callsR=c(NULL),callsF=c("outlier","outlierp"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mcorr",infile="misc.r",visible=FALSE,callsR=c("corrrisk"),callsF=c("mcorr"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dti3Dreg",infile="misc.r",visible=FALSE,callsR=c(NULL),callsF=c("dti3Dreg"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dti3Dev",infile="misc.r",visible=FALSE,callsR=c(NULL),callsF=c("dti3Dev"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dti3Dall",infile="misc.r",visible=FALSE,callsR=c(NULL),callsF=c("dti3Dall"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dti3Dand",infile="misc.r",visible=FALSE,callsR=c(NULL),callsF=c("dti3Dand"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dtieigen",infile="misc.r",visible=FALSE,callsR=c(NULL),callsF=c("dtieigen"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dtiind3D",infile="misc.r",visible=FALSE,callsR=c(NULL),callsF=c("dtiind3D"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="kldist",infile="misc.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="fwhm2bw",infile="misc.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="replind",infile="misc.r",callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="sofmchi",infile="misc.r",callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="fncchir",infile="misc.r",callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="fncchis",infile="misc.r",callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="fncchiv",infile="misc.r",callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="fncchiL",infile="misc.r",callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="fncchiL2",infile="misc.r",callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="Spatialvar.gauss",infile="misc.r",callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="Varcor.gauss",infile="misc.r",callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="corrrisk",infile="misc.r",visible=TRUE,callsR=c("thcorr3D"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="thcorr3D",infile="misc.r",visible=TRUE,callsR=c(NULL),callsF=c("thcorr"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="andir2.image",infile="misc.r",visible=TRUE,callsR=c(NULL),callsF=c("thcorr"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="andir.image",infile="misc.r",visible=TRUE,callsR=c(NULL),callsF=c("thcorr"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="connect.mask",infile="misc.r",visible=TRUE,callsR=c(NULL),callsF=c("lconnect"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="sphcoord",infile="misc.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="gettriangles",infile="misc.r",visible=FALSE,callsR=c(NULL),callsF=c("distvert","triedges"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="create.designmatrix.dti",infile="misc.r",callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="identifyFA",infile="misc.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mthomogen",infile="misc.r",visible=FALSE,callsR=c("extract"),callsF=c("mthomog","parofor"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="vcrossp",infile="misc.r",visible=FALSE,callsR=c("extract"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="showFAColorScale",infile="misc.r",visible=FALSE,callsR=c("extract"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
#
#    File mixtensorpl.r
#
list(funct="dwiMixtensor",infile="mixtensorpl.r",visible=TRUE,callsR=c("dtieigen","sioutlier","getsiind3","getsiind2","selisample","plmatrix","pgetsiind2","pmixtens","pmixtns0","pmixtns1","pmixtns2"),callsF=c("sweeps0","sweeps0p"),callsC=c("mixture","mixtrl0","mixtrl1","mixtrl2"),purpose="",moveto=""), 
list(funct="selisample",infile="mixtensorpl.r",visible=TRUE,callsR=c(NULL),callsF=c("selisamp"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getsiind3",infile="mixtensorpl.r",visible=TRUE,callsR=c("selisample","plmatrix","pgetsii30","pgetsii31"),callsF=c("iandir","getsii30","getsii31"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getsiind2",infile="mixtensorpl.r",visible=TRUE,callsR=c("selisample"),callsF=c("getsii"),callsC=c(NULL),purpose="",moveto=""), 
#
#    File improve.r
#
list(funct="mfunpl",infile="improve.r",visible=TRUE,callsR=c(NULL),callsF=c("mfunpl"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="gmfunpl",infile="improve.r",visible=TRUE,callsR=c(NULL),callsF=c("gmfunpl"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mfunwghts",infile="improve.r",visible=TRUE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mfunwghtsi",infile="improve.r",visible=TRUE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mfunpl0",infile="improve.r",visible=TRUE,callsR=c(NULL),callsF=c("mfunpl0"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mfunpl0h",infile="improve.r",visible=TRUE,callsR=c(NULL),callsF=c("mfunpl0h"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="gmfunpl0",infile="improve.r",visible=TRUE,callsR=c(NULL),callsF=c("mfunpl0g"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mfunplwghts0",infile="improve.r",visible=TRUE,callsR=c(NULL),callsF=c("mfunpl0","mfunpl0w"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mfunplwghts0h",infile="improve.r",visible=TRUE,callsR=c(NULL),callsF=c("mfunpl0h"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dwiMtImprove",infile="improve.r",visible=TRUE,callsR=c("extract","mfunpl0","gmfunpl0","mfunpl0h","mfunpl","gmfunpl","mfunpli","gmfunpli","mfunplwghts0","mfunplwghts0h","mfunwghts","mfunwghtsi"),callsF=c("sweepimp","imprparb"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dwiMtCombine",infile="improve.r",visible=TRUE,callsR=c("extract"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getcube3",infile="improve.r",visible=TRUE,callsR=c(NULL),callsF=c("mfunpl"),callsC=c(NULL),purpose="",moveto=""), 
#
#    File sqrtODF.r
#
list(funct="dwiSqrtODF",infile="sqrtODF.r",visible=TRUE,callsR=c("sioutlier","connect.mask","Lnlm","getD0","Kofgrad","mcorr"),callsF=c("sqrteap","Mofcall"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="Qlm",infile="sqrtODF.r",visible=TRUE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="kappan",infile="sqrtODF.r",visible=TRUE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="hnni",infile="sqrtODF.r",visible=TRUE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="g1f1",infile="sqrtODF.r",visible=TRUE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="Inna",infile="sqrtODF.r",visible=TRUE,callsR=c("hnni","g1f1","kappan"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="Yab",infile="sqrtODF.r",visible=TRUE,callsR=c("sphcoord","getsphericalharmonicseven"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="Qlmmat",infile="sqrtODF.r",visible=TRUE,callsR=c("Qlm"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="QlmmatR",infile="sqrtODF.r",visible=TRUE,callsR=c("Qlmmat"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="Kofgrad",infile="sqrtODF.r",visible=TRUE,callsR=c("Yab","QlmmatR","Inna"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="Lnlm",infile="sqrtODF.r",visible=TRUE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="Mofcall",infile="sqrtODF.r",visible=TRUE,callsR=c("Lnlm"),callsF=c("Mofcall"),callsC=c(NULL),purpose="",moveto=""), 
#
#    File aws.r
#
list(funct="dti.smooth",infile="aws.r",visible=TRUE,callsR=c("dtilin.smooth","dtireg.smooth"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dtilin.smooth",infile="aws.r",visible=FALSE,callsR=c("sioutlier","dtiTensor","Spatialvar.gauss","andir2.image","dti3Dev"),callsF=c("projdt2","smsigma","awssidti"),callsC=c(NULL),purpose="",moveto=""), 
#
#    File aws2.r
#
list(funct="dtireg.smooth",infile="aws2.r",visible=FALSE,callsR=c("sdpar","sioutlier","dtiTensor","andir2.image","getvofh","gethani","Spatialvar.gauss"),callsF=c("projdt2","awsrgdti"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="gethani",infile="aws2.r",visible=FALSE,callsR=c(NULL),callsF=c("gethani"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getvofh",infile="aws2.r",visible=FALSE,callsR=c(NULL),callsF=c("getvofh"),callsC=c(NULL),purpose="",moveto=""), 
#
#    File ricebias.r
#
list(funct="dwiRiceBias",infile="ricebias.r",visible=TRUE,callsR=c("ricebiascorr"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="ricebiascorr",infile="ricebias.r",visible=FALSE,callsR=c("sofmchi"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
#
#    File fibers.r
#
list(funct="selectFibers",infile="fibers.r",visible=TRUE,callsR=c(NULL),callsF=c("roifiber"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="reduceFibers",infile="fibers.r",visible=TRUE,callsR=c("sortFibers"),callsF=c("reducefe","reducefi"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="ident.fibers",infile="fibers.r",visible=FALSE,callsR=c("sortFibers"),callsF=c("fibersta"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="sortFibers",infile="fibers.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="compactFibers",infile="fibers.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="expandFibers",infile="fibers.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
#
#    File tracking.r
#
list(funct="tracking",infile="tracking.r",visible=TRUE,callsR=c("dtiIndices","ident.fibers","compactFibers","extract"),callsF=c(NULL),callsC=c("interface_tracking","interface_tracking_mixtensor"),purpose="",moveto=""), 
#
#    File parallel.r
#
list(funct="pmatrix",infile="parallel.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="plmatrix",infile="parallel.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="pnlrdtirg",infile="parallel.r",visible=FALSE,callsR=c(NULL),callsF=c("nlrdtirp"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="pnltens",infile="parallel.r",visible=FALSE,callsR=c("opttensR","rho2D","tensRres"),callsF=c("nlrdtirp"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="pmixtens",infile="parallel.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c("mixture"),purpose="",moveto=""), 
list(funct="pmixtns2",infile="parallel.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c("mixtrl2"),purpose="",moveto=""), 
list(funct="pmixtns1",infile="parallel.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c("mixtrl1"),purpose="",moveto=""), 
list(funct="pmixtns0",infile="parallel.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c("mixtrl0"),purpose="",moveto=""), 
list(funct="pgetsii30",infile="parallel.r",visible=FALSE,callsR=c(NULL),callsF=c("pgtsii30"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="pgetsii31",infile="parallel.r",visible=FALSE,callsR=c(NULL),callsF=c("pgtsii31"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="pgetsiind2",infile="parallel.r",visible=FALSE,callsR=c(NULL),callsF=c("getsii"),callsC=c(NULL),purpose="",moveto=""), 
# 
#   File varest.r
#
list(funct="awssigmc",infile="varest.r",visible=TRUE,callsR=c("sofmchi","fncchis"),callsF=c("paramw3","awsvchi","awsadchi"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="afsigmc",infile="varest.r",visible=TRUE,callsR=c(NULL),callsF=c("afmodevn","afmodem1"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="awslinsd",infile="varest.r",visible=TRUE,callsR=c("setawsdefaults","awsgfamily","gethani","Spatialvar.gauss","awsgsigma2"),callsF=c("caws03d","cgaws"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="awsgfamily",infile="varest.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="awsgsigma2",infile="varest.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="setawsdefaults",infile="varest.r",visible=FALSE,callsR=c("getvofh"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="IQRdiff",infile="varest.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
# 
#   File awsse3.r
#
list(funct="dwi.smooth",infile="awsse3.r",visible=TRUE,callsR=c("getmask","awssigmc","sofmchi","getnext3g","suggestkappa","getkappasmsh","gethseqfullse3msh","getkappas","gethseqfullse3","interpolatesphere","lkfullse3msh","fncchiv","reduceparam","create.designmatrix.dti"),callsF=c("adsmse3m","adsmse3p","asmse30p"),callsC=c(NULL),purpose="Note: actual Multishell implementation is still method dwi.smooth.ms",moveto=""), 
list(funct="lkfullse3",infile="awsse3.r",visible=FALSE,callsR=c(NULL),callsF=c("lkfulse3"),callsC=c(NULL),purpose="",moveto=""), 
# 
#   File smse3misc.r
#
list(funct="betagamma",infile="smse3misc.r",visible=FALSE,callsR=c(NULL),callsF=c("bgstats"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="matrm",infile="smse3misc.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getkappas",infile="smse3misc.r",visible=FALSE,callsR=c("betagamma"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="suggestkappa",infile="smse3misc.r",visible=FALSE,callsR=c("getkappas"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="vredsphere",infile="smse3misc.r",visible=FALSE,callsR=c("getkappas"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
# 
#   File paramse3.r
#
list(funct="gethseqfullse3",infile="paramse3.r",visible=FALSE,callsR=c(NULL),callsF=c("ghfse3i"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="reduceparam",infile="paramse3.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
# 
#   File multishell.r
#
list(funct="sphtrarea1",infile="multishell.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getsphwghts",infile="multishell.r",visible=FALSE,callsR=c("sphtrarea1"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getnext3g",infile="multishell.r",visible=FALSE,callsR=c("getsphwghts"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="interpolatesphere",infile="multishell.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="needs to be replaced by interpolatesphere0",moveto=""), 
list(funct="lkfullse3msh",infile="multishell.r",visible=FALSE,callsR=c("lkfullse3"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="gethseqfullse3msh",infile="multishell.r",visible=FALSE,callsR=c("gethseqfullse3"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getkappasmsh",infile="multishell.r",visible=FALSE,callsR=c("getkappas"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
# 
#   File awsse3ms.r  
#
list(funct="dwi.smooth.ms",infile="awsse3ms.r",visible=TRUE,callsR=c("getmask","awssigmc","sofmchi","getnext3g","getkappasmsh3","gethseqfullse3msh","interpolatesphere0","lkfullse3msh","lkfulls0","create.designmatrix.dti"),callsF=c("adsmse3c"),callsC=c(NULL),purpose="Note: actual Multishell implementation ",moveto=""), 
list(funct="getkappas3",infile="awsse3ms.r",visible=FALSE,callsR=c("betagamma"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getkappasmsh3",infile="awsse3ms.r",visible=FALSE,callsR=c("getkappas3"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="interpolatesphere0",infile="multishell.r",visible=FALSE,callsR=c(NULL),callsF=c("ipolsp"),callsC=c(NULL),purpose="",moveto=""), 
list(funct="lkfulls0",infile="multishell.r",visible=FALSE,callsR=c(NULL),callsF=c("lkfuls0"),callsC=c(NULL),purpose="",moveto=""), 
# 
#   File dki.r  
#
list(funct="dkiTensor",infile="dki.r",visible=TRUE,callsR=c("dkiDesign","connect.mask","pseudoinverseSVD"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dkiDesign",infile="dki.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dkiIndices",infile="dki.r",visible=TRUE,callsR=c("dkiDesign","rotateKurtosis","kurtosisFunctionF1","kurtosisFunctionF2"),callsF=c("dti3DevAll"),callsC=c(NULL),purpose="",moveto=""), 
# 
#   File dkimisk.r  
#
list(funct="rotateKurtosis",infile="dkimisc.r",visible=FALSE,callsR=c("defineKurtosisTensor"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="defineKurtosisTensor",infile="dkimisk.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="kurtosisFunctionF1",infile="dkimisc.r",visible=FALSE,callsR=c("kurtosisFunctionF2"),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="kurtosisFunctionF2",infile="dkimisc.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="pseudoinverseSVD",infile="dkimisc.r",visible=FALSE,callsR=c(NULL),callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""),
#
#  File awsse3mstestprop.r
#
list(funct="dwi.smooth.testprop",infile="awsse3mstestprop.r",visible=TRUE,c("awssigmc","sofmchi","getnext3g","getkappasmsh3","gethseqfullse3msh","interpolatesphere0","lkfullse3msh","lkfulls0","create.designmatrix.dti")
,callsF=c("adsmse3c","exceed"),callsC=c(NULL),purpose="",moveto=""))

#
#   Fortran subroutines
#

Fstructure <- list(
list(funct="iandir",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getsii30",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getsii31",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getsii",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="sweeps0",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="sweeps0p",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="initdata",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="nlrdtirg",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="opttensR",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="tensRres",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="rho2D0",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="D2rho0",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dti2Dga",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dti2Dfa",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getmask",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="outlier",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="outlierp",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mcorr",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dti3Dreg",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dti3Dall",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dti3Dev",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dti3Dand",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dtieigen",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="dtiind3D",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="thcorr",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="lconnect",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="distvert",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="triedges",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mthomog",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="parofor",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mfunpl",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="gmfunpl",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mfunpl0",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mfunpl0h",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mfunpl0g",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mfunpli",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="gmfunpli",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mfunplwghts0",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mfunplwghts0h",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mfunwghts",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mfunwghtsi",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="smsigma",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="awssidti",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="projdt2",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="awsrgdti",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="roifiber",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="reducefe",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="reducefi",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="fibersta",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="pgtsii31",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="pgtsii30",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="gethani",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="getvofh",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="paramw3",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="awsvchi",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="awsadchi",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="afmodevn",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="afmodem1",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="caws03d",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="cgaws",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="ghfse3i",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="bgstats",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="adsmse3m",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="adsmse3p",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="asmse30p",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="lkfulse3",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="ipolsp",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="lkfuls0",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""),
list(funct="exceed",infile="misc.f",callsF=c(NULL),callsC=c(NULL),purpose="",moveto="")))


#
#   C/C++ routines
#

Cstructure <- list(
list(funct="mixture",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mixtrl0",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mixtrl1",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="mixtrl2",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="interface_tracking",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="interface_tracking_mixtensor",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""), 
list(funct="",infile="",callsF=c(NULL),callsC=c(NULL),purpose="",moveto=""))

callsofR <- function(funct,struct=Rstructure){
n <- length(struct)
flist <- NULL
for(i in 1:n){
if(funct %in% struct[[i]]$callsR){
flist <-c(flist,struct[[i]]$funct)
}
}
flist
}
callsofRall <- function(struct=Rstructure){
n <- length(struct)
flist <- NULL
for(i in 1:n){
flist <- c(flist,paste(struct[[i]]$funct,"called by ",paste(callsofR(struct[[i]]$funct),collapse=", ")))
}
flist
}
callsofRall()
