pdf(file='dmso_test.pdf',width=4.5,height=4.5);
gstable=read.table('dmso_test.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("TNFAIP3","PRDM1","HAT1","RNF168","CHD3","SUZ12","UBE2B","HDAC2","ALKBH3","H3F3B")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='3,4,5_vs_2 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2267.5288660825545,1773.110483579203,1443.0737428748512,1833.9927853004212),c(1165.1153851867123,440.3866798889651,637.5728621013569,1396.875333881668),c(1438.7236948195996,1133.2488742875776,811.0369693011472,864.7323495457945),c(2683.641503649237,890.4098297973824,1602.6976287555094,3513.093951945936),c(3586.548925437766,3388.18285883939,3721.1741720040122,1903.361424329919),c(949.6488413508134,474.11432495704776,400.4437368335585,964.5091591087707),c(2705.3021614951745,2166.2784603728524,2016.0589054869247,2496.3207497327494),c(1193.6162507734714,823.9181866631623,1360.032414960058,854.2795409249112),c(1875.356955608749,1843.456714721204,2562.286306882009,2158.0298525478015),c(1650.7701347850873,1227.686280478209,1085.9960328412403,1683.8524432913712))
targetgene="TNFAIP3"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(510.73551131472317,198.51128240071512,514.8562330717181,304.0817053347849),c(2200.266823297803,1623.745198277694,1798.3060900659111,1774.1266995626354),c(2836.406143194266,2443.8087969330754,1777.0844173765752,2830.810625601013),c(2040.661976011952,1882.0025947990127,1766.0122403212695,1997.436701917868),c(421.8128106840348,254.40280851353782,161.46924872320906,725.9950714867989),c(4026.602290097326,3561.6393191895295,3367.787187655503,2804.2034763842194),c(1705.4917967116648,1412.706504851691,1259.4601400410306,1103.2464371677663),c(2083.9832917038257,1403.070034832239,1086.9187142625158,2101.9647881267006),c(3381.3426932131,3019.106057094371,2208.8993225335,2718.680496758811),c(4581.799151727394,3865.1881248022737,2891.6835742773555,4313.208939108089))
targetgene="PRDM1"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1002.0904340304502,1363.560507752485,1475.367592619493,1394.0245678941544),c(3034.7721676781093,2089.186700217235,1865.6618338190212,1738.967252383301),c(2489.8356176592756,943.4104149043694,1042.6300060412927,751.6519653744214),c(3757.5541189583205,2057.386349153043,2356.528349937577,2999.0058188643156),c(1776.1739433668274,2365.753389775513,1991.1465071124867,2739.5861140005773),c(4234.088591568932,3747.6231905649565,3297.6633996385667,2657.8641556918537),c(3678.8917299388654,2132.5508153047695,1887.8061879296329,3060.7724152604437),c(1358.921271176674,1387.6516828011154,1476.2902740407685,1249.5857578601315),c(3705.1125262786836,3284.1089826293064,3947.231120216505,2601.799091270753),c(2501.2359638939793,1694.091429419695,1667.285328244793,2217.895938285587))
targetgene="HAT1"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2339.3510473611873,766.0993665464491,1843.5174797084098,1324.6559288646565),c(3054.1527562771057,2676.047724401873,2093.564144874065,2332.876833115303),c(1441.0037640665403,2124.8416392892077,2264.260207810029,2043.999213047257),c(2524.0366563633866,2188.4423414175926,2064.961020814525,2217.895938285587),c(3838.4965772247165,4132.118344341099,2433.110907903442,3119.6882456690587),c(796.8842018057846,881.7370067798754,1236.3931045091438,1261.939077139357),c(899.4873179181174,502.0600880134591,339.5467630293768,151.09059733822122),c(1504.8457029808808,572.4063191554601,205.75795694443212,342.09191850163296),c(320.34972919517236,352.69480271195016,420.74272810161904,398.15698292273396),c(2429.413782615346,2253.006690547922,1313.8983438962841,1295.1980136603493))
targetgene="RNF168"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2465.894890566398,1724.928133481942,1889.6515507721838,2394.643429511431),c(965.6093260793984,1220.9407514645925,1160.7332279645543,998.718350958934),c(2766.864031162574,2513.191381073131,3373.3232761831564,2408.8972594489987),c(1944.8990676404414,1538.9442621065148,1583.3213189087244,2496.3207497327494),c(6396.734272292214,1409.8155638458554,1989.3011442699355,1273.3421410894116),c(2787.3846543850404,1520.6349690695556,3257.0654171024457,2535.281218228769),c(92.34280450109951,308.3670406224701,60.896973804181705,448.52051536880765),c(884.6668678130027,2181.696812403976,1021.4083333519568,1238.1826939100772),c(2154.6654383589885,2722.3027804952435,1174.5734492836866,2453.5592599200454),c(54.72166192657748,340.1673916866623,86.73205359989515,0.0))
targetgene="CHD3"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1695.2314851004317,1668.0729603671741,1338.810742270722,2233.1000235523265),c(2679.081365155356,2310.825510664635,1699.5791779894348,1661.0463153912624),c(1565.2675380248102,2008.2403520538364,1054.6248645178741,817.2195830872344),c(2016.7212489190742,2744.4666615399838,2821.559786260419,2126.671426685152),c(1157.1351428224198,914.5010048460128,443.80976363350607,1112.7489904594784),c(2780.544446644218,1983.1855300032607,2512.461510133133,2370.8870462821505),c(1898.1576480781564,1908.9847108534789,559.1449412929411,1242.9339705559332),c(1475.2048027706512,1099.521229219495,1430.1562029769946,1192.5704381098594),c(4606.879913443742,3944.207178961781,4372.587255424502,3147.2456502150235),c(4776.7450723408265,4967.600295027604,4163.138572794967,3539.70110116273))
targetgene="SUZ12"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2754.3236503044,1923.4394158826572,1637.7595227639777,2425.0516000449093),c(2139.8449882538735,1402.1063878302937,1743.8678862106578,2945.7915204307283),c(3074.673379499572,2242.406573526525,2637.023502005323,2888.776200680456),c(3347.141654508989,2664.4839603785304,2765.2762195626146,2013.5910425137786),c(2339.3510473611873,1981.2582359993703,2666.5493074861383,1332.257971498026),c(2543.4172449623825,1787.5651886083813,2997.7919377240355,2422.2008340573957),c(2050.922287623185,2203.860693448716,1615.615168653366,1568.8715484616557),c(1037.4315073580315,1055.1934671300148,1193.9497591304716,1605.9315062993326),c(2504.65606776439,4045.3901141660294,1613.7698058108151,1946.1229141426231),c(4105.264679116782,3497.0749700591996,2104.6363219293708,2960.9956056974675))
targetgene="UBE2B"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2153.525403735518,2442.8451499311304,1730.9503463128012,1526.1100586489515),c(2443.0941980969906,2347.4440967385535,1782.6205059042281,1951.8244461176503),c(3267.3392308660636,2943.9415909426443,2416.502642320483,2839.3629235635535),c(3578.5686830734735,2629.7926683085025,3168.4880006599997,2702.5261561629004),c(4334.411638434325,3618.4944923042976,2837.245370422102,4269.497193966214),c(3138.5153184139126,2606.665140261817,2117.5538618272276,2669.2672196419085),c(1387.4221367634332,980.992647980233,1539.0326106875013,1849.1968705671604),c(4001.5215283809785,2650.9929023512977,2492.1625188650723,3197.6091826610973),c(1055.6720613335572,2489.100206024501,1982.8423743210074,1584.075633728395),c(5095.954766912528,2707.8480754660654,4739.814461092143,3860.8874024225966))
targetgene="HDAC2"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1919.8183059240932,2265.53410157321,2030.8218082273322,2130.4724480018363),c(2403.192986275528,1731.6736624955586,1038.939280356191,1932.8193395342264),c(855.0259676027731,861.5004197390258,529.6191358121257,1184.9683954764898),c(3165.8761493772013,2127.7325802950436,2605.6523336819564,3236.5696511571164),c(1071.6325460621424,361.36762572945713,476.1036133781479,572.053708161064),c(3988.9811475228043,1912.8392988612597,1706.9606293596387,2554.286324812193),c(2184.306338569218,1288.3960416007578,1312.052981053733,1694.3052519122543),c(1746.533043156598,1106.2667582331114,1994.8372327975885,1729.464699091589),c(2840.966281688148,1703.7278994391472,1594.39349596403,2683.5210495794763),c(4196.467448994411,5232.60322056254,4952.953869406779,5191.24486326228))
targetgene="ALKBH3"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2009.8810411782522,1604.4722582387897,1567.6357347470412,1422.5322277692903),c(2633.479980216541,1870.43883077567,1403.3984417600057,1589.777165703422),c(84.36256213680696,186.9475183773725,256.5054351145835,365.848301730913),c(1399.9625176216073,904.8645348265607,847.9442261521665,1269.541119772727),c(2083.9832917038257,2344.5531557327176,1478.1356368833196,3136.79284159414),c(2123.8845035252884,1663.254725357448,2084.3373306613103,1784.5795081835188),c(1948.3191715108524,1670.0002543710646,1332.3519723217937,1355.064099398135),c(1832.0356399168752,1530.2714390890078,2074.18783502728,1225.8293746308516),c(2728.1028539645818,1432.9430918925407,1587.9347260151017,1714.2606138248498),c(1450.1240410543032,1169.867460361496,1033.403191828538,704.1391989158612))
targetgene="H3F3B"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("RRM1","CHD9","RAN","SMN1","SHPRH","EIF4A3","SETD4","METTL14","PSMA3","ARCN1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='3,4,5_vs_2 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(0.0,41.43682108364442,0.0,0.0),c(326.0499023125242,37.58223307586354,35.984575429743735,189.10081050506935),c(0.0,0.0,46.134071063774016,0.0),c(84.36256213680696,6.745529013616533,106.10836344668024,111.17987351303071),c(206.3462668481359,145.51069729372807,403.2117810973849,131.13523542562598),c(0.0,0.0,130.09808039984273,0.9502553291712027),c(0.0,84.80093617117927,0.0,62.71685172529938),c(0.0,24.091175048630475,85.80937217861967,305.03196066395606),c(0.0,0.0,0.0,161.54340595910446),c(0.0,44.327762089480075,0.0,85.52297962540824))
targetgene="RRM1"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(462.8540571289679,132.98328626844022,313.71168323366334,114.03063950054432),c(0.0,245.72998549603085,0.0,208.10591708849338),c(5183.737432919746,3642.585667352928,4429.793503543581,3698.393741134321),c(0.0,601.3157292138167,790.7379780330866,321.1863012598665),c(1030.5912996172092,1237.3227504976612,1630.3780713937738,980.6634997046812),c(1922.0983751710341,2239.515632520689,3688.8803222593706,2565.6893887622473),c(3521.5669518999553,2408.1538578611026,2758.817449613686,4463.349281117139),c(1495.7254259931178,2862.995242779246,3578.1585517063127,1235.3319279225634),c(2547.977383456264,2069.9137601783304,1411.702574551485,2224.5477255897854),c(1219.8370471132896,2198.0788114370444,1257.6147771984797,1842.545083262962))
targetgene="CHD9"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(664.6401854832224,304.5124526146892,199.29918699550376,385.8036636435083),c(0.0,0.963647001945219,0.0,0.0),c(63.841938914340396,62.637055126439236,356.1550286123354,0.0),c(220.0266823297803,302.58515861079877,119.02590334453697,141.5880440465092),c(1.1400346234703642,142.61975628789241,223.28890394866625,186.25004451755572),c(624.7389736617596,1319.2327456630048,465.9541177441176,966.4096697671132),c(1537.9067070615213,1674.8184893807907,1633.1461156576001,1411.129163819236),c(0.0,88.65552417896015,23.98971695316249,112.13012884220191),c(278.16844812676885,53.964232108932265,454.8819406888118,100.72706489214748),c(1491.1652874992365,1406.9246228400198,699.3925173268141,867.583115533308))
targetgene="RAN"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(361.39097564010547,931.8466508810268,1060.1609530455269,360.1467697558858),c(69.54211203169221,220.67516344545515,262.0415236422364,130.18498009645478),c(947.3687721038726,450.0231499084173,638.4955435226324,338.29089718494816),c(26.220796339818378,482.7871479745547,471.4902062717705,16.154340595910448),c(308.9493829604687,280.42127756605873,431.8149051569248,394.3559616060491),c(1961.999586992497,2694.3570174388324,3469.282143995806,2132.372958660179),c(241.68734017571722,971.3561779607808,1239.1611487729701,668.0294964073555),c(345.43049091152034,225.49339845518125,0.0,147.28957602153642),c(2195.7066848039217,3366.01897779465,3668.58133099131,2043.0489577180858),c(425.23291455444587,360.4039787275119,723.3822342799766,725.9950714867989))
targetgene="SMN1"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(534.6762384076009,2014.985881067453,710.4646943821199,974.961967729654),c(755.8429553608514,1223.8316924704282,1086.9187142625158,302.18119467644243),c(1446.703937183892,1244.0682795112778,1300.0581225771518,1220.1278426558242),c(901.767387165058,465.4415019395408,1082.3053071561385,1173.5653315264353),c(876.6866254487101,1234.4318094918256,1839.8267540233078,1421.5819724401192),c(57.00173117351821,210.07504642405775,270.34565643371576,423.8138768103564),c(672.6204278475149,956.9014729316025,489.0211532760046,962.6086484504283),c(1023.751091876387,2527.6460861023097,1604.5429915980603,1437.7363130360297),c(2896.8279782381956,4110.918110298304,2851.085591741234,3873.240721701822),c(1931.2186521587969,2336.843979717156,3130.658062387705,2838.4126682343826))
targetgene="SHPRH"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2806.765242984037,2319.4983336821424,1777.0844173765752,2541.9330055329674),c(79.80242364292549,57.81882011671314,56.283566697804304,180.5485125425285),c(2.2800692469407284,114.67399323148106,65.5103809105591,363.9477910725706),c(1049.9718882162053,419.1864458461703,162.39193014448455,579.6557507944336),c(1.1400346234703642,102.14658220619322,0.9226814212754804,96.92604357546267),c(63.841938914340396,102.14658220619322,176.23215146361676,234.71306630528707),c(327.1899369359945,403.76809381504677,130.09808039984273,114.03063950054432),c(933.6883566222283,782.4813655795178,1229.01165313894,1311.3523542562598),c(79.80242364292549,0.0,58.128929540355266,30.408170533478486),c(0.0,0.0,54.43820385525334,0.0))
targetgene="EIF4A3"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1269.9985705459858,2783.9761886197375,2637.023502005323,1126.052565067875),c(1037.4315073580315,639.8616092916254,990.0371650285904,536.8942609817295),c(1055.6720613335572,789.2268945931344,865.4751731564006,995.8675849714205),c(1252.8980511939303,1583.2720241959948,3371.477913340605,1701.907294545624),c(859.5861060966546,1316.3418046571692,1828.754576968002,1809.28614674197),c(1492.3053221227067,1959.0943549546303,2899.065025647559,2426.9521107032515),c(2837.5461778177364,5069.746877233797,3222.0035230939775,3026.5632234102804),c(1000.9503994069797,1363.560507752485,2105.559003350646,2463.0618132117575),c(897.2072486711767,1640.1271973107628,2627.7966877925683,1284.745205039466),c(2529.7368294807384,2883.2318298200953,4090.2467405142042,2237.8513001981823))
targetgene="SETD4"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(368.23118338092763,517.4784400445826,290.6446477017763,690.8356243074644),c(920.007941140584,589.751965190474,1406.166486023832,462.7743453063757),c(1474.064768147181,1993.7856470246581,3969.3754743271165,3520.695994579306),c(637.2793545199336,1754.8011905422438,1669.130691087344,1501.4034200905003),c(2555.9576258205566,2985.3784120262885,2525.3790500309897,4928.024137081858),c(397.8720835911571,176.34740135597508,729.8410042289049,1104.1966924969374),c(1300.7795053796856,2512.227734071186,2447.8738106438495,1700.0067838872817),c(2801.065069866685,3166.5440483919897,1878.579373716878,4279.950002587097),c(853.8859329793028,878.8460657740397,2131.3940831463597,1357.9148653856487),c(1103.5535155193127,1615.072375260187,2069.5744279209025,1561.2695058282861))
targetgene="METTL14"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1022.6110572529167,300.65786460690833,272.19101927626673,568.2526868443792),c(1584.6481266238063,1054.2298201280696,1077.691900049761,907.4938393584986),c(905.1874910354692,680.3347833733246,1169.0373607560336,1475.7465262028777),c(176.70536663790645,284.2758655738396,544.3820385525335,809.6175404538646),c(1244.9178088296378,1452.216031931445,894.0782972159404,1837.793806617106),c(379.6315296156313,184.05657737153683,419.8200466803436,453.2717920146637),c(539.2363769014822,407.62268182282764,641.2635877864589,200.50387445512376),c(5.700173117351821,296.80327659912746,7.381451370203843,228.06127900108865),c(1615.4290614575061,1000.2655880191373,1974.538241529528,723.1443054992852),c(0.0,22.163881044740037,60.896973804181705,62.71685172529938))
targetgene="PSMA3"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(525.5559614198379,335.3491566769362,386.60351551442625,300.28068401810003),c(511.8755459381935,577.2245541651862,813.8050135649737,496.983537156539),c(497.0550958330788,633.1160802780089,429.04686089309837,707.940220232546),c(3.4201038704110927,157.0744613170707,107.9537262892312,55.11480909192976),c(0.0,43.364115087534856,362.61379856126376,452.3215366854925),c(86.64263138374768,99.25564120035756,52.59284101270238,375.35085502262507),c(197.225989860373,158.03810831901592,201.14454983805473,374.40059969345384),c(349.9906294054018,276.56668955827786,231.59303674014558,554.9491122359824),c(175.56533201443608,0.963647001945219,0.0,56.065064421100956),c(324.9098676890538,44.327762089480075,97.80423065520093,101.6773202213187))
targetgene="ARCN1"
collabel=c("H2171.DMSO.18.R2.0","H2171.DMSO.18.R3.0","H2171.PRMT5.18.R1.0","H2171.PRMT5.18.R2.0")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("dmso_test_summary.Rnw");
library(tools);

texi2dvi("dmso_test_summary.tex",pdf=TRUE);

