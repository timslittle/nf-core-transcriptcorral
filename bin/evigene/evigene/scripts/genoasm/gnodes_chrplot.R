#
# chromosome plots for gnodes DNA coverage depth, xCopy deviations and feature densities
#
# source("gnodes_chrplot.R")
#
# data prep, see evigene/scripts/genoasm/gnodes_covsum2p.pl
#   aweed20gnodes/
#   pt=arath18tair_chr; 
#   pt=arath20max_chr; 
#   env kucg=33 $evigene/scripts/genoasm/gnodes_covsum2p.pl -plotchr -debug \
#    -genexcopy arath18tair1cds_SRR10178325_b2_mim.geneycopy -asmid $pt  -title${pt}_testplota -anntab $pt.anntab  \
#    -crclass arath18tair1cds.idclass -sumdata arath20asm.metad   ${pt}_SRR10178325_test8f.covtab
#   
#   env kucg=0 $evigene/scripts/genoasm/gnodes_covsum2p.pl -plotchr -debug  \
#    -genexcopy arath18tair1cds_SRR3703081_b2_mim.geneycopy  -asmid $pt  -title ${pt}_testplhetz -anntab $pt.anntab  \
#    -crclass arath18tair1cds.idclass -sumdata arath20asm.metad   ${pt}_SRR3703081_hetest8f.covtab
#   
#   dromel20gnodes/
#   pt=drosmel6ref_chr;
#   pt=drosmel20pi_chr;  
#   env kucg=0 $evigene/scripts/genoasm/gnodes_covsum2p.pl -plotchr -debug \
#     -genexcopy dromel6relt1cds_SRR11460802_b2_mim.genexcopy -asmid $pt  -title ${pt}_testplota \
#    -anntab $pt.anntab  -crclass dromel6relt1cds.idclass -sumdata drosmelchr.metad ${pt}_SRR11460802_b2_mim.cc8a.covtab
# 
#   env kucg=0 $evigene/scripts/genoasm/gnodes_covsum2p.pl -plotchr -debug \
#     -genexcopy dromel6relt1cds_SRR10512945_1odd.genexcopy   -asmid $pt  -title ${pt}_testplhetz\
#    -anntab $pt.anntab  -crclass dromel6relt1cds.idclass -sumdata drosmelchr.metad ${pt}_SRR10512945_1odd_bwa.cc8a.covtab
#----------------

try2dmel = function() {

  dtag="drosmel20pi8a";
  dfn="crtabs/drosmel20pi_testplota_SRR11460802.plottab"; 
  gnodes_chrplot(dtag,dfn);

  dtag="drosmel6ref8a";
  dfn="crtabs/drosmel6ref_testplota_SRR11460802.plottab"; 
  gnodes_chrplot(dtag,dfn);

  dtag="drosmel20hetz8b";
  dfn="crtabs/drosmel20pi_testplhetz_SRR10512945.plottab"; 
  gnodes_chrplot(dtag,dfn);
  
  dtag="dropse20uc8b";
  dfn="crtabs/dropse20chrs_testplota_SRR11813283a.plottab";
  gnodes_chrplot(dtag,dfn);
}
 
try2arath = function() {

  dtag="arath18tair8b";
  dfn="crtabs/arath18tair8_testplota_SRR10178325.plottab";
  gnodes_chrplot(dtag,dfn);

  dtag="arath20max8b";
  dfn="crtabs/arath20max8_testplota_SRR10178325.plottab";
  gnodes_chrplot(dtag,dfn);

  dtag="arath18bwa8b";
  dfn="crtabs/arath18tair8_testplbwa_SRR10178325.plottab";
  gnodes_chrplot(dtag,dfn);

  dtag="arath18hetz8b";
  dfn="crtabs/arath18tair8_testplhetz_SRR3703081.plottab";
  gnodes_chrplot(dtag,dfn);
}

#UPD21JAN: xcopytrans="n" default, ="l" or "log" for log transform, change xcopymax for that?
gnodes_chrplot =  function(dtag,dfn,pvers="xcopy", #chrxcopy
  crmin=100, binmax=500, crscale=1e6, xcopymax=3, xcopytrans="n", eoswap=F, pstyle="s", yrev=T, chrorder="size" ) {
  
  ctab = read.table(dfn,header=T,stringsAsFactors=F)
  
  ccols= grep("ChrID|ObsLoc|EstLoc|Xall|Nall|Nuniq|Ndup|NCDSann|NTEann|NRPTann|NNOann|NCONTAM", colnames(ctab));
  pltab=ctab[,ccols];
  # daph: count CONTAM, remove all-contam scafs  
  icon=grep("NCONTAM",colnames(ctab)); iall=grep("Nall",colnames(ctab)); # which name? Nallasm or Nall
  if( icon>0 & iall>0 ) { 
    rcon= ctab[,icon] >= ctab[,iall]; 
    if(sum(rcon)>0) { pltab= ctab[ ! rcon, ccols ]; }
  }
  if(eoswap) { pltab[,2:3]= pltab[,3:2]; dtag=paste(dtag,"Est",sep="_"); }
  plcolnames= colnames(pltab)
  plcolnames=gsub("asm","",gsub("ann","",plcolnames));
  
  jcols= grep("Xall|Ndup|NCDS|NTE|NRPT", plcolnames)
  # check jcols have vals .. drop if all zero; jcolnames == plotcolnames should be option
  # FIXME? Nall has small vals for small scafs wedged into chrasm.  *Should* normalize Nxxx by Nall,ie
  #  for j in jcols { if(plcolnames[j] =~ /^N/) { pltab[j]= round(100*pltab[j]/ pltab[iall]); } }

  # OPTION for log(Xall) transform to show *high* extremes : 
  #  if(Xall>1) Xall=1+ln(Xall); == 10>3.3, 20>4, 50>5, 100>5.6, 200>6.3
  if(xcopytrans == "l") {
     ix=grep("Xall",colnames(pltab));
     jx=(pltab[,ix] > 1); pltab[jx,ix]= 1 + log(pltab[jx,ix]);
     #? reset xcopymax? not this: max(pltab[jx,ix]), maybe 90% quartile? but want same xcopymax over sets, leave
     # jx=(pltab[,ix] < 1); # no change for low end?
  }
  
  colmax=rep(binmax,ncol(pltab)); names(colmax)=plcolnames;
  colclr=rainbow(length(jcols),start=0.3)
  names(colclr)=plcolnames[jcols]; 
  colclr["Ndup"]= "orange"; colclr["Xall"]= "red";
  
  clabel=plcolnames; clabel=gsub("asm","",gsub("ann","",clabel)); 
  clabel=gsub("^N","",clabel,perl=T);
  clabel=gsub("dup","DUP",clabel);
  clabel=gsub("Xall","XCOPY",clabel);
  if(xcopytrans == "l") { clabel=gsub("XCOPY","lXCOPY",clabel); } #?
  
  cbaseline=rep(NA,ncol(pltab)); names(cbaseline)=plcolnames; 
  colmax["Xall"] = xcopymax; # true max=6 for dmel6ref, lower for dmel20pi
  cbaseline["Xall"]= 1;
  plotname= paste(dtag,pvers,sep="_")
  
  # BUG in chrplotrot: wont do > 10 plots? .. 
  # FIXME for >10, up to 20..30 cr > crmin, do several plots, 
  # order crlist by crn sizes, not name ie not Chr1 Chr10 .. Chr19 Chr2 ..
  crn= table(pltab[,1]); 
  cror= order( crn, decreasing=T ); 
  if(chrorder == "size") { crn= crn[cror]; }
  # if(chrorder == "name") { crors= order( names(crn), decreasing=F ); crn= crn[crors]; }
  crok= crn >= crmin; 
  #xx if(sum(crok) < 2 & length(crn) >= 2) { ci=min(length(crn),9); crmin= crn[cror][ci]; crok= crn >= crmin; }
  crlist= names(crn)[crok];
  if(chrorder == "name") {  crlist= sort(crlist ); } 
  
  maxcrs=50; # 50? debug at 9
  crlista= crlist;
  if(length(crlista)>maxcrs){ crlista= crlista[ 1:maxcrs]; } # max plots=5, 10 chr each
  while( length(crlista) > 0) {
    crlist= crlista; 
    if(length(crlista)>10){ crlist= crlista[ 1:10]; crlista= crlista[ -c(1:10) ]; }
    else { crlista=c(); }
    
    cmaxb=0; for (chr in crlist) { cmaxb= max( cmaxb, pltab[ (pltab[,1] == chr), 2] ); }
  
    icr=0; ncr= length(crlist); 
    for (chr in crlist) { 
     icr= 1+icr; tip= ncr*10 + icr;
     # new plotname.chr for each icr==1
     chrplotrot(plotname, pltab, chr, icols=jcols, triplot=tip, idcol=0,
        colbaseline=cbaseline, chrmaxbp=cmaxb, colnam=clabel, colclr=colclr, colmax=colmax, pstyle=pstyle, yrev=yrev, crscale=crscale);
    }
  }
  
#----  
#   if(length(crlist)>10) { crlist= crlist[1:10]; }
#   
#   cmaxb=0; for (chr in crlist) { cmaxb= max( cmaxb, pltab[ (pltab[,1] == chr), 2] ); }
# 
#   icr=0; ncr= length(crlist); 
#   for (chr in crlist) {
#    icr= 1+icr; tip= ncr*10 + icr;
#    chrplotrot(plotname, pltab, chr, icols=jcols, triplot=tip, idcol=0,
#       colbaseline=cbaseline, chrmaxbp=cmaxb, colnam=clabel, colclr=colclr, colmax=colmax, pstyle=pstyle, yrev=yrev);
#   # was: title=paste(dtag,chr)  : dtag too long, need only chr default
#   }
  
}


chrplotrot = function(plnam, pltab, chr, icols, sppnam="",vers="8g",
   triplot=0, title=NA, colnam=NA, colmax=NA, colclr=NA, colbaseline=NA, colisdensity=NA,
   idcol=0, chrmaxbp = 0, pwid=0, pstyle='l', smoo=F, yrev=F, crscale=1e6) {

  bandht= 16 
  ncols= length(icols)
  cset <- (pltab[,1] == chr)  
  chrlab=chr;
  if(length(grep("^\\w+:\\w",chr,perl=T))>0) { chrlab=gsub("^\\w+:","",chr,perl=T); } # drop prefix: asmid:Chr 

  if(is.na(title)) { title= paste(sppnam,chrlab); }
  
  #NOTE: x is now vert, y is horiz
  ymax <- bandht * ncols
  xmin <- min(pltab[cset,2]) # for subset spans, eg 5 Mb .. 10 Mb
  xmax <- max(pltab[cset,2]) # expect 10,000,000 to 15,000,000 for 10 daph chrs
  if(chrmaxbp < 100) { chrmaxbp=xmax; } else { xmax = max(chrmaxbp, xmax); }
  
  if(pwid == 0) pwid <- max(8, 1 * round((xmax-xmin)/1000000) ) # 1-inch per 1 Mbp, for 10-20 Mbp chrs
  ylim= c(xmin,xmax); if(yrev) ylim= c(xmax,xmin);
  
  # FIXME: triplot calc bug for nplot > 10, change opts: triplot > iplot=1..n, nplot=n, 
  if(triplot >= 20) {
    nplot= trunc( ifelse( triplot>=100, 10 + (triplot - 100)/20, triplot/10 ) ); 
    plen = 3*nplot; 
    triploti= triplot - nplot*10; triplote= triploti >= nplot
  } else {
    plen <- ifelse(triplot > 0,  9, min(8, ncols) )
    nplot= ifelse(triplot>0, 3, 1); 
    triploti= triplot; triplote= triploti >= nplot
  }
  
  # if(anyNA(colclr)) { colclr <- rainbow(ncols,start=0.3); }
  clriset <- colclr; 
  if(anyNA(colnam)) { colnam=colnames(pltab) }

  if(smoo) { for (i in icols) { pltab[cset,i]= runmed(pltab[cset,i], 3); } }
  # for (i in icols) { if(is.na(colmax[i])) colmax[i]= max(pltab[cset,i]); }
  
  if(anyNA(colisdensity)) {
  colisdensity= rep( FALSE, ncol(pltab)) # test new "s" plot, not for xcopy/colbaseline 
  for (ic in icols) { if( is.na(colbaseline[ic]) ) { colisdensity[ic]= TRUE; } }
  }
  
  ret=triplot
  pltfile=paste(plnam,vers,chrlab,".pdf",sep="")
  if(triplot < 2 | triploti == 1) {
    ret=pltfile
    pdf( pltfile,height=pwid,width=plen) 
    par( mar=c(3.1,3.1,3.1,3.1), cex=0.85)  #was cex=0.75
    if(triplot > 0) { par( mfcol=c(1,nplot) ) }
  }
  
  plot( pltab[cset,icols[1]], pltab[cset,2], 
    type="n", xaxt="n", yaxt="n", xlim=c(0,ymax), ylim=ylim,  frame.plot=F,
    main="", sub = "", xlab="", ylab="")
  title( main=title, sub = title, line=0, cex.main = 1.2, cex.sub = 1.2) # was line=1 

  xmaxc <- max(pltab[cset,2]) ;  xminc <- min(pltab[cset,2]) ; 
  xtic1 <- seq(trunc(xminc/crscale)*crscale,trunc(xmaxc/crscale)*crscale,crscale) #? xmin
  if(crscale < 1e6) {
  xwidc= (xmaxc - xminc)/1000000;
  if(length(xtic1)/xwidc > 20) {
   crs=crscale*10;
   xtic1 <- seq(trunc(xminc/crs)*crs,trunc(xmaxc/crs)*crs,crs) #? xmin
  }
  }
  xtics <- axTicks(2) 
  xtics <- xtics[xtics <= xmaxc]
  lxM= ifelse(crscale < 2e4, "K","M"); # fixme C,K,M 10K? 100K?
  xtlab <- paste(round(xtics/crscale),lxM) ; # fixme M for diff crscale != 1e6
    if(yrev) { xtlab[length(xtlab)]= "1"; } else { xtlab[1] <- "1"; }
  axis(2, at=xtic1, labels=F, tcl=-0.3)
  if(length(xtics) == length(xtlab)) {  # fixme bug for diff crscale..
  axis(2, at=xtics, labels=xtlab)
  }

  if(idcol > 0) {
    cids= cset & (pltab[,idcol] != '.')  # fixme '.' == na
    if(sum(cids,na.rm=T) > 0) {
    xv <- pltab[cids,2]
    # remove overlap ids here? if(xv[i] - xv[i-1] < minx) cids[i]=F; 
    idv <- pltab[cids,idcol]
    text( 1, xv, labels=idv, cex=0.60)
    }
  }
  
  ci<-1; for (ic in icols) {  
    coly0 <- 4 + ci*bandht - bandht # steps at 20,40,60,80,100 = top
    xv <- pltab[cset,2]
    yv <- pltab[cset,ic]
    yv <- ( bandht/colmax[ic]) * yv; # need colmin ?
    yv[ yv > bandht ] = bandht; yv[ yv < 0 ] = 0;     
    np <- length(yv)

    if( ! is.na(colbaseline[ic]) ) { 
      yb= bandht * colbaseline[ic]/colmax[ic] ; 
      points( x=c(coly0+yb,coly0+yb), y=c(xminc,xmaxc), type="l", col = "black"); 
    }
    
    if(pstyle == 'p') { points( coly0+yv, xv, col=clriset[ci], pch = "."); }
    if(pstyle == 'l') { points( coly0+yv, xv, col=clriset[ci], type="l" ); }
    #o: if(pstyle == 's') { segments( rep(coly0,np), xv,  coly0+yv, xv, col=clriset[ci]); }
    # this new s is good, for density, percent graph cols; need to ensure colisdensity false for others
    if(pstyle == 's') { 
      if(colisdensity[ic]) { segments( rep(coly0,np), xv,  coly0+yv, xv, col=clriset[ci]); }
      points( coly0+yv, xv, col=clriset[ci], type="l" );
    }
    
    p1=1; p3=3; if(yrev) { p1=3; p3=1; }
    text( coly0+4, xminc, colnam[ic], pos=p1) # , offset = 0.1, cex=1.1 
    text( coly0+4, xmaxc, colnam[ic], pos=p3) # , offset = 0.1, cex=1.1 
    ci<- ci+1
  } 
   
  if(triplot == 0 || triplote) { dev.off() }
  return(ret);
}

