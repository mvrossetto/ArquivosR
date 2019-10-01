print("Loading libraries")
suppressWarnings(suppressMessages(library(stats4)))
suppressWarnings(suppressMessages(library(SummarizedExperiment)))
suppressWarnings(suppressMessages(library(TCGAbiolinks)))
suppressWarnings(suppressMessages(library(survival)))
suppressWarnings(suppressMessages(library(survcomp)))
suppressWarnings(suppressMessages(library(scales)))
suppressWarnings(suppressMessages(library(dunn.test)))
suppressWarnings(suppressMessages(library(agricolae)))
suppressWarnings(suppressMessages(library(gdata)))
print("Loaded libraries")
#library(dplyr)
#library(DT)

args = commandArgs(trailingOnly=TRUE)

args[1] = 'D:\\ArquivosR\\'
args[2] = 'TCGA-STAD'
args[3] = 10
args[4] = 1
args[5] = "OLFML2B"
#args[6] = ""
# args[7] = "TERT"
# args[8] = "TFF1"
# args[9] = "TFF2"

caminho = args[1]     #'D:/ArquivosR/' #pegar de argumento
NameProject = args[2]
QtdPackpageForGet = as.numeric(args[3])
CaminhoCompleto = paste(caminho,NameProject,sep="")
#print(CaminhoCompleto)

#Salva parametros recebidos
write(args,file = paste(CaminhoCompleto,"/parametros.txt",sep = ""),sep="\n")

## Testa se o caminho existe, se sim, carrega i .rdata já com os dados baixados
## Se não, precisa fazer todo o download e salvar o .rdata para a próxima vez.
if (file.exists(paste(CaminhoCompleto,"data.RData",sep="/"))){
  print("Loading data")
  isfar<-load(paste(CaminhoCompleto,"data.RData",sep="/")) 
  setwd(CaminhoCompleto)
  
  samplesNT <- TCGAquery_SampleTypes(getResults(query,cols="cases"),"NT")
  samplesTP <- TCGAquery_SampleTypes(getResults(query,cols="cases"),"TP")
  
  
  if ((length(samplesNT) <= 0) || (length(samplesNT)) <= 0)
    stop("Don't have Normal data")
  
  
  print("Loaded data")
} else {
  setwd(caminho)
  
  query <- GDCquery(project = NameProject, #parametro
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    platform = "Illumina HiSeq", 
                    file.type  = "results",
                    experimental.strategy = "RNA-Seq",
                    legacy = TRUE) 
  
  
  print("Getting data of TCGA")
  TCGAbiolinks::GDCdownload(query, directory = caminho,
                            files.per.chunk = QtdPackpageForGet, method = "api")  
  stad.exp <- GDCprepare(query, directory = caminho)
  stad.assay <- assay(stad.exp, "raw_count")
  
  
  stad.assay.prep <- TCGAanalyze_Preprocessing(stad.exp, cor.cut = 0.6)
  
  stad.assay.norm <- TCGAanalyze_Normalization(tabDF = stad.assay.prep,
                                               geneInfo =  geneInfo,
                                               method = "gcContent")
  
  stad.assay.filt <- TCGAanalyze_Filtering(tabDF = stad.assay.norm,
                                           method = "quantile", 
                                           qnt.cut =  0.25) 
  
  samplesNT <- TCGAquery_SampleTypes(getResults(query,cols="cases"),"NT")
  samplesTP <- TCGAquery_SampleTypes(getResults(query,cols="cases"),"TP")
  
  
  if ((length(samplesNT) <= 0) || (length(samplesNT)) <= 0)
    stop("Don't have Normal data")
  
  
  #  dataDEGs.edgeR <- TCGAanalyze_DEA(mat1 = stad.assay.filt[,samplesNT],
  #                                    mat2 = stad.assay.filt[,samplesTP],
  #                                  Cond1type = "Non-tumor",
  #                                  Cond2type = "Tumor",
  #                                 fdr.cut = 1,
  #                                 logFC.cut = 0,
  #                                 pipeline = "edgeR",
  #                                 method = "glmLRT")
  
  save.image(paste(CaminhoCompleto,"data.RData",sep="/"))
  #dev.off()
}

#Joga novamente os valore nas variaveis para ter certeza que as mesmas estão com os valores corretos
#Pois as mesmas podem estar salvas no .rdata carregado a cima.
args = read.table(paste(CaminhoCompleto,"/parametros.txt",sep = ""), header = F)
args = as.character(args[,1])

caminho = args[1]     #'D:/ArquivosR/' #pegar de argumento
NameProject = args[2]   
QtdPackpageForGet = as.numeric(args[3])
CaminhoCompleto = paste(caminho,NameProject,sep="")

QtdGenes = as.numeric(args[4])
QtdGenes =  4 + as.numeric(QtdGenes)

#print(QtdGenes)
NameOfGenes = c()
q = 5
for (q in 5: QtdGenes){
  NameOfGenes <-  append(NameOfGenes,args[q])
}


### PCA --------
names <- colnames(stad.assay.filt)
names(names) <- substr(names, 14, 16)
color.stad <- ifelse(names(names)=="01A","brown3", "dodgerblue3")
stad.assay.pc <- stad.assay.filt
stad.assay.pc <- data.matrix(stad.assay.pc)
transp <- t(stad.assay.pc)
## analyse principal component
p3A <- prcomp(transp, retx=TRUE, center=TRUE, scale=TRUE)
pca.stad <- p3A$x

NameOfImage = "PCA TCGA-STAD.jpg"

DirectoryOfImage = paste(CaminhoCompleto,NameOfImage,sep="/")

png(file=DirectoryOfImage,
    width=800, height=600)  
# tiff(file=DirectoryOfImage,res = 300)


plot(pca.stad[,1], pca.stad[,2], xlab="STAD PCA 1", ylab="STAD PC 2",
     type="p", cex.lab=1.2, cex.axis=1.2, 
     xlim=c(-80,150), ylim=c(-160,170),
     col=color.stad,
     pch=16) #main="STAD", cex.main=1.2, font.main=1, 
legend("bottomleft", legend=c("non-GC (n=35)", "GC (n=415)"), bty="n",
       cex=1.2,col=c('dodgerblue3','brown3'), y.intersp=1.2, pch=16)

dev.off()

### END PCA



## STAD
# boxplot(stad.norm1, las=2, cex.lab=1, cex.axis=0.9, col=group.stad$color, #par(font.main=1),
#         outpch=20, outcol=group.stad$color, outcex=0.9, #main="STAD", cex.main=1.4,
#         ylab='STAD gene expression', xaxt="n", par(font.lab=1))

#Repeticao para cada Genes
print("Starting data analysis")
i = 1
for (i in 1:length(NameOfGenes)) { 
  
  codicao = NameOfGenes[i] %in% row.names(stad.assay.filt)
  if (codicao){
    
    Genes_stad_edgeR <- as.data.frame(stad.assay.filt[NameOfGenes[i],])
    colnames(Genes_stad_edgeR) <-  "Genes.assay.filt"
    
    Genes_stad_edgeR$Genes.assay.filt = Genes_stad_edgeR$Genes.assay.filt + 1 
    
    Genes_stad_edgeR$genes.log2.assay.filt <- log2(Genes_stad_edgeR$Genes.assay.filt)
    
    Group = factor(stad.exp$shortLetterCode,levels = c("NT","TP"))
    Genes_stad_edgeR$Group <- Group  
    
    #head(Genes_stad_edgeR$Group)
    testTP = shapiro.test(subset(Genes_stad_edgeR$genes.log2.assay.filt, Genes_stad_edgeR$Group == "TP"))
    testNT = shapiro.test(subset(Genes_stad_edgeR$genes.log2.assay.filt, Genes_stad_edgeR$Group == "NT"))
    
    if ((testNT$p.value < 0.05 ) | (testTP$p.value < 0.05)){
      Teste = wilcox.test(Genes_stad_edgeR$genes.log2.assay.filt ~ Genes_stad_edgeR$Group )
    }else{
      Teste = t.test(Genes_stad_edgeR$genes.log2.assay.filt ~ Genes_stad_edgeR$Group )
    }
    
    mean.stad = tapply(Genes_stad_edgeR$genes.log2.assay.filt,Genes_stad_edgeR$Group,mean)
    logfc = mean.stad["TP"] -mean.stad["NT"]   
    
    
    #Troca Titulo boxplot
    Label = NameOfGenes[i]
  
    if (Teste$p.value < 0.001) {
      #Não precisa fazer nada.
      #Teste$p.value = format(Teste$p.value,scientific = TRUE)
    }else{
      Teste$p.value = round(Teste$p.value,3)
    }
    
    NameOfImage = paste(NameOfGenes[i],".jpg",sep="")
    
    
    
    if (!file.exists(paste(CaminhoCompleto,NameOfGenes[i],sep="/"))){
      dir.create(file.path(paste(CaminhoCompleto,NameOfGenes[i],sep="/")))
    }
    
    DirectoryOfImage = paste(paste(CaminhoCompleto,NameOfGenes[i],sep="/"),NameOfImage,sep="/")
    
    png(file=DirectoryOfImage,
       width=800, height=600)

    #tiff(file=DirectoryOfImage,res = 300)
    par(mar = c(7,5,3,2))
    bx_color<-c("cornflowerblue","brown3")
    boxplot(genes.log2.assay.filt ~ Group, 
            data=Genes_stad_edgeR, varwidth=F,
            ###main= bquote(expression(italic(.(Label)))), troque pra por ali embaixo no title###
            medlwd=2,
            plot=T,outline=T,outcol="white",par(font.lab=1),
            names=c("Adjacent non-tumoral tissue","Gastric cancer"),
            las=1,
            boxwex=0.40,par(cex.lab=1.3),par(cex.axis=1.3),
            boxcol=bx_color, whiskcol=bx_color,
            medcol=bx_color,staplecol=bx_color,
            col=alpha(bx_color,0.5))
    
    title(main=bquote(italic(.(NameOfGenes[i]))), cex.main=1.3, font.main=1)
    
    axis(side=1, at=1:2, labels=c("",""),lty=1,las=1)
    title(ylab = "Relative gene expression (log2)", cex.lab = 1.3, line = 3.5)
    
    tumor= Genes_stad_edgeR[Genes_stad_edgeR$Group == 'TP' ,]
    tumor= tumor[tumor$genes.log2.assay.filt != '-Inf' ,]
    tumor= na.omit(tumor)
    normal = Genes_stad_edgeR[Genes_stad_edgeR$Group == 'NT' ,]
    normal= normal[normal$genes.log2.assay.filt != '-Inf' ,]
    normal= na.omit(normal)
    
    #Monta a legenda inferior
    legendBottomNT <- paste("median = ",round(median(normal$genes.log2.assay.filt),2)   ,
                            "\n mean = ",round(mean(normal$genes.log2.assay.filt),2) ,
                            "\n n = ",length(samplesNT) )
    legendBottomTP <- paste("median = ",round(median(tumor$genes.log2.assay.filt),2), 
                            "\n mean = ",round(mean(tumor$genes.log2.assay.filt),2) ,
                            "\n n = ",length(samplesTP))
    
    mtext(c(legendBottomNT,legendBottomTP),
          side=1,line=5, at=1:2, cex=1.3)
    
    stripchart(genes.log2.assay.filt ~ Group,
               vertical=TRUE, data=Genes_stad_edgeR, 
               method="jitter", add=TRUE, pch=19, col=alpha(bx_color,0.8))
    
    legendTopLogFC = "";
    if (logfc < 0.001){
      legendTopLogFC <- paste('LogFC = ', format(logfc,scientific = T))
    }else{
      legendTopLogFC <- paste('LogFC = ', round(logfc,3))
    }
    
    
    legendTop = "";
    #valorLegendTop = "";
    if (Teste$p.value < 0.001){
      legendTop <- paste(legendTop,'\nP = ', format(Teste$p.value,scientific = T))
    }else{
      legendTop <- paste(legendTop,'\nP = ', round(Teste$p.value,3))
    }
    
    legend("topleft", legendTopLogFC, bty="n",cex=1.3)
    legend("topleft", legendTop, bty="n",cex=1.3)
    
    
    dev.off()
    
    
    # Start clinical analyze --------------------------------------------------
    stad.clinical_test <- stad.exp@colData
    stad.clinical = as.data.frame(stad.clinical_test)
    
    stad.clinical <- cbind(stad.clinical, Genes_stad_edgeR$genes.log2.assay.filt)
    
    Genes.log2.subtypes <- Genes_stad_edgeR
    Genes.log2.subtypes <- cbind(Genes.log2.subtypes, stad.clinical$barcode)
    colnames(Genes.log2.subtypes)[length(Genes.log2.subtypes)] = "barcode"
    Genes.log2.subtypes <- cbind(Genes.log2.subtypes, stad.clinical$subtype_Lauren.Class)
    colnames(Genes.log2.subtypes)[length(Genes.log2.subtypes)] = "subtype_Lauren.Class"
    Genes.log2.subtypes <- cbind(Genes.log2.subtypes, stad.clinical$subtype_Molecular.Subtype)
    colnames(Genes.log2.subtypes)[length(Genes.log2.subtypes)] = "subtype_Molecular.Subtype"
    Genes.log2.subtypes <- cbind(Genes.log2.subtypes, stad.clinical$subtype_WHO.Class)
    colnames(Genes.log2.subtypes)[length(Genes.log2.subtypes)] = "subtype_WHO.Class"
    Genes.log2.subtypes <- cbind(Genes.log2.subtypes, stad.clinical$subtype_KRAS.mutation)
    colnames(Genes.log2.subtypes)[length(Genes.log2.subtypes)] = "subtype_KRAS.mutation"
    Genes.log2.subtypes <- cbind(Genes.log2.subtypes, stad.clinical$subtype_TP53.mutation)
    colnames(Genes.log2.subtypes)[length(Genes.log2.subtypes)] = "subtype_TP53.mutation"
    Genes.log2.subtypes <- cbind(Genes.log2.subtypes, stad.clinical$subtype_PIK3CA.mutation)
    colnames(Genes.log2.subtypes)[length(Genes.log2.subtypes)] = "subtype_PIK3CA.mutation"
    Genes.log2.subtypes <- cbind(Genes.log2.subtypes, stad.clinical$subtype_ARID1A.mutation)
    colnames(Genes.log2.subtypes)[length(Genes.log2.subtypes)] = "subtype_ARID1A.mutation"
    Genes.log2.subtypes <- cbind(Genes.log2.subtypes, stad.clinical$subtype_RHOA.mutation)
    colnames(Genes.log2.subtypes)[length(Genes.log2.subtypes)] = "subtype_RHOA.mutation"
    Genes.log2.subtypes <- cbind(Genes.log2.subtypes, stad.clinical$subtype_ARHGAP26.ARHGAP6.CLDN18.Rearrangement)
    colnames(Genes.log2.subtypes)[length(Genes.log2.subtypes)] = "subtype_ARHGAP26.ARHGAP6.CLDN18.Rearrangement"
    Genes.log2.subtypes <- cbind(Genes.log2.subtypes, stad.clinical$gender)
    colnames(Genes.log2.subtypes)[length(Genes.log2.subtypes)] = "gender"
    Genes.log2.subtypes <- cbind(Genes.log2.subtypes, stad.clinical$days_to_death)
    colnames(Genes.log2.subtypes)[length(Genes.log2.subtypes)] = "days_to_death"
    Genes.log2.subtypes <- cbind(Genes.log2.subtypes, stad.clinical$days_to_last_follow_up)
    colnames(Genes.log2.subtypes)[length(Genes.log2.subtypes)] = "days_to_last_follow_up"
    
    
    ##Tumor stage filter
    Genes.log2.subtypes <- cbind(Genes.log2.subtypes, stad.clinical$tumor_stage)
    colnames(Genes.log2.subtypes)[length(Genes.log2.subtypes)] = "tumor_stage"
    Genes.log2.subtypes$tumor_stage = gsub("stage ","",tumorstage$tumor_stage)
    Genes.log2.subtypes$tumor_stage = gsub("a","",tumorstage$tumor_stage)
    Genes.log2.subtypes$tumor_stage = gsub("b","",tumorstage$tumor_stage)
    Genes.log2.subtypes$tumor_stage = gsub("c","",tumorstage$tumor_stage)
    
    
    
    Genes.log2.lauren <- Genes.log2.subtypes
    
    Genes.log2.lauren1 <- subset(Genes.log2.lauren, Genes.log2.lauren$subtype_Lauren.Class =="Diffuse")
    Genes.log2.lauren1 <- rbind(Genes.log2.lauren1, subset(Genes.log2.lauren, Genes.log2.lauren$subtype_Lauren.Class=="Intestinal"))
    Genes.log2.lauren1 <- rbind(Genes.log2.lauren1, subset(Genes.log2.lauren, Genes.log2.lauren$subtype_Lauren.Class=="Mixed"))
    group_NT <- subset(Genes.log2.lauren, Genes.log2.lauren$Group=="NT")
    group_NT$subtype_Lauren.Class<-"NT"
    group_NT$subtype_Molecular.Subtype<-"NT"
    group_NT$subtype_WHO.Class<-"NT"
    Genes.log2.lauren1 <- rbind(Genes.log2.lauren1, subset(group_NT, group_NT$subtype_Lauren.Class=="NT"))
    
    Genes.log2.lauren1$subtype_Lauren.Class = as.character(Genes.log2.lauren1$subtype_Lauren.Class)
    Genes.log2.lauren1$subtype_Lauren.Class = as.factor(Genes.log2.lauren1$subtype_Lauren.Class)
    
    # Normality test
    STDif = shapiro.test(subset(Genes.log2.lauren1$genes.log2.assay.filt, Genes.log2.lauren1$subtype_Lauren.Class == "Diffuse"))
    STInt = shapiro.test(subset(Genes.log2.lauren1$genes.log2.assay.filt, Genes.log2.lauren1$subtype_Lauren.Class == "Intestinal"))
    STMix = shapiro.test(subset(Genes.log2.lauren1$genes.log2.assay.filt, Genes.log2.lauren1$subtype_Lauren.Class == "Mixed"))
    STnt  = shapiro.test(subset(Genes.log2.lauren1$genes.log2.assay.filt, Genes.log2.lauren1$subtype_Lauren.Class == "NT"))
    
    Diffuse = Genes.log2.lauren1[Genes.log2.lauren1$subtype_Lauren.Class == 'Diffuse',]
    Intestinal = Genes.log2.lauren1[Genes.log2.lauren1$subtype_Lauren.Class == 'Intestinal',]
    Mixed = Genes.log2.lauren1[Genes.log2.lauren1$subtype_Lauren.Class == 'Mixed',]
    NT = Genes.log2.lauren1[Genes.log2.lauren1$subtype_Lauren.Class == 'NT',]
    
    qtdDiffuse = nrow(Diffuse)
    qtdIntestinal = nrow(Intestinal)
    qtdMixed = nrow(Mixed)
    qtdNt = nrow(NT)
    
    #Se for não paramétrico = mediana se não média
    if ((STDif$p.value < 0.05 ) | (STInt$p.value < 0.05)| (STMix$p.value < 0.05)| (STnt$p.value < 0.05)){
      ## Non parametric Kruskal.Wallis Test
      kw.stad.lauren <- kruskal.test(genes.log2.assay.filt ~ subtype_Lauren.Class, data=Genes.log2.lauren1)
      
      GroupLetrinhas <- kruskal(Genes.log2.lauren1$genes.log2.assay.filt,
                                Genes.log2.lauren1$subtype_Lauren.Class,
                                alpha=0.025, group=TRUE, p.adj="BH")$groups
      
      mdif = round(median(Diffuse$genes.log2.assay.filt),2)
      mint = round(median(Intestinal$genes.log2.assay.filt),2)
      mmix = round(median(Mixed$genes.log2.assay.filt),2)
      mnt = round(median(NT$genes.log2.assay.filt),2)
      
      legendaDif = paste("Median =",mdif,"\nn =",qtdDiffuse)
      legendaInt = paste("Median =",mint,"\nn =",qtdIntestinal)
      legendaMix = paste("Median =",mmix,"\nn =",qtdMixed)
      legendaNt = paste("Median =",mnt,"\nn =",qtdNt)  
      
      
      if (kw.stad.lauren$p.value < 0.001){
        legtopright =  paste('P = ', format(kw.stad.lauren$p.value,scientific = T))
      }else{
        legtopright = paste("P = ",(round(kw.stad.lauren$p.value,3)))
      }
      #Não paramétrico
      TipoTeste = "Kruskal-Wallis with Benjamini-Hochberg correction"
    }else{ 
      anTULP3<-aov(genes.log2.assay.filt ~ subtype_Lauren.Class, data=Genes.log2.lauren1)
      
      GroupLetrinhas <- HSD.test(anTULP3, 'subtype_Lauren.Class', alpha = 0.05, group=TRUE)$groups
      
      mdif = round(mean(Diffuse$genes.log2.assay.filt),2)
      mint = round(mean(Intestinal$genes.log2.assay.filt),2)
      mmix = round(mean(Mixed$genes.log2.assay.filt),2)
      mnt = round(mean(NT$genes.log2.assay.filt),2)
      
      legendaDif = paste("Mean =",mdif,"\nn =",qtdDiffuse)
      legendaInt = paste("Mean =",mint,"\nn =",qtdIntestinal)
      legendaMix = paste("Mean =",mmix,"\nn =",qtdMixed)
      legendaNt = paste("Mean =",mnt,"\nn =",qtdNt)  
      
      
      p.anova =  as.data.frame(unlist(summary(anTULP3)))
      
      if (p.anova['Pr(>F)1',] < 0.001){
        legtopright =  paste('P = ', format(p.anova['Pr(>F)1',],scientific = T))
      }else{
        legtopright = paste("P = ",(round(p.anova['Pr(>F)1',],3)))
      }
    
      TipoTeste = "Anova with Benjamini-Hochberg correction"
    }
    
    NameOfImage = paste(NameOfGenes[i]," expression among Lauren Classification.jpg",sep="")
    
    DirectoryOfImage = paste(paste(CaminhoCompleto,NameOfGenes[i],sep="/"),NameOfImage,sep="/")
    
    png(file=DirectoryOfImage,
        width=800, height=600)  
    #tiff(file=DirectoryOfImage,res = 300)
    bx_color4 <- c("#FF6E00","brown3","mediumorchid3","cornflowerblue")
    
    Label = NameOfGenes[i]
    par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1 ## bottom, left, top, and right
    boxplot(genes.log2.assay.filt ~ subtype_Lauren.Class, 
            data=Genes.log2.lauren1, varwidth=F,
            #main=bquote(expression(italic(.(Label)))), Colocado abaixo
            medlwd=2,
            plot=T,outline=T,outcol="White",par(font.lab=1),
            names=c("Diffuse","Intestinal","Mixed","Non-tumor"),
            #ylim=c(9.3,13.3),
            las=1,
            boxwex=0.40,par(cex.lab=1.3),par(cex.axis=1.3),
            boxcol=bx_color4, whiskcol=bx_color4,
            medcol=bx_color4,staplecol=bx_color4,
            col=alpha(bx_color4,0.5))
    
    axis(side=1, at=1:4, labels=c("","","",""),lty=1,las=1)
    
    title(main=bquote(italic(.(NameOfGenes[i]))), cex.main=1.3, font.main=1)
    title(ylab = "Relative gene expression (log2)", cex.lab = 1.3, line = 3.5)
    
    mtext(c(legendaDif,legendaInt,legendaMix,legendaNt),side=1,line=3.5, at=1:4, cex=1.3)
    
    
    stripchart(genes.log2.assay.filt ~ subtype_Lauren.Class,
               vertical=TRUE, data=Genes.log2.lauren1, 
               method="jitter", add=TRUE, pch=19, col=alpha(bx_color4,0.8))
    
    
    legend("topright",legend=bquote(italic(.(legtopright))), bty="n",cex=1.3)
    
    
    #letrinhas  "Diffuse","Intestinal","Mixed","Non-tumor"
    #### Como verificar quantos nives tem
    if ( table(GroupLetrinhas[,2])[1] != 4 ){
      text(x=1, y=max(Diffuse$genes.log2.assay.filt), GroupLetrinhas["Diffuse",2],pos=3, cex=1.3)
      text(x=2, y=max(Intestinal$genes.log2.assay.filt), GroupLetrinhas["Intestinal",2],pos=3, cex=1.3)
      text(x=3, y=max(Mixed$genes.log2.assay.filt), GroupLetrinhas["Mixed",2],pos=3, cex=1.3)
      text(x=4, y=max(NT$genes.log2.assay.filt), GroupLetrinhas["NT",2],pos=3, cex=1.3)
    }
    
    dev.off()
    
    
    
    ## Tumor Stage Comparing ----
    Genes.log2.TumorStage = Genes.log2.subtypes
    Genes.log2.TumorStage <- subset(Genes.log2.TumorStage, Genes.log2.TumorStage$tumor_stage!="not reported")  
    
    stageI = Genes.log2.TumorStage[Genes.log2.TumorStage$tumor_stage == 'i',]
    stageIi = Genes.log2.TumorStage[Genes.log2.TumorStage$tumor_stage == 'ii',]
    stageIii = Genes.log2.TumorStage[Genes.log2.TumorStage$tumor_stage == 'iii',]
    stageIv = Genes.log2.TumorStage[Genes.log2.TumorStage$tumor_stage == 'iv',]
    
    qtdstageI = nrow(stageI)
    qtdstageIi = nrow(stageIi)
    qtdstageIii = nrow(stageIii)
    qtdstageIv = nrow(stageIv)
    
    ## shapiro
    STstageI = shapiro.test(subset(Genes.log2.TumorStage$genes.log2.assay.filt, Genes.log2.TumorStage$tumor_stage == "i"))
    STstageIi = shapiro.test(subset(Genes.log2.TumorStage$genes.log2.assay.filt, Genes.log2.TumorStage$tumor_stage == "ii"))
    STstageIii  = shapiro.test(subset(Genes.log2.TumorStage$genes.log2.assay.filt, Genes.log2.TumorStage$tumor_stage == "iii"))
    STstageIv = shapiro.test(subset(Genes.log2.TumorStage$genes.log2.assay.filt, Genes.log2.TumorStage$tumor_stage == "iv"))
    
    Genes.log2.TumorStage$tumor_stage = as.factor(Genes.log2.TumorStage$tumor_stage)
    
    if ((STstageI$p.value < 0.05 ) | (STstageIi$p.value < 0.05)| (STstageIii$p.value < 0.05)| (STstageIv$p.value < 0.05)){
      ## Non parametric Kruskal.Wallis Test
      kw.stad.tumor_stage <- kruskal.test(genes.log2.assay.filt ~ tumor_stage, data=Genes.log2.TumorStage)
      
      mstageI = round(median(stageI$genes.log2.assay.filt),2)
      mstageIi = round(median(stageIi$genes.log2.assay.filt),2)
      mstageIii = round(median(stageIii$genes.log2.assay.filt),2)
      mstageIv = round(median(stageIv$genes.log2.assay.filt),2)
      
      
      legendastageI = paste("Median =",mstageI,"\nn =",qtdstageI)
      legendastageIi = paste("Median =",mstageIi,"\nn =",qtdstageIi)
      legendastageIii = paste("Median =",mstageIii,"\nn =",qtdstageIii)
      legendastageIv = paste("Median =",mstageIv,"\nn =",qtdstageIv)
      
      
      
      if (kw.stad.tumor_stage$p.value < 0.001){
        legtoprightTumorStage =  paste('P = ', format(kw.stad.tumor_stage$p.value,scientific = T))
      }else{
        legtoprightTumorStage = paste("P = ",(round(kw.stad.tumor_stage$p.value,3)))
      }
      
      GroupLetrinhasTumorStage <- kruskal(Genes.log2.TumorStage$genes.log2.assay.filt,
                                         Genes.log2.TumorStage$tumor_stage,
                                         alpha=0.025, group=TRUE, p.adj="BH")$groups
      
      #Não paramétrico
      TipoTeste = "Kruskal-Wallis with Benjamini-Hochberg correction"
    }else{ 
      ### Parametric t test
      antumorstage<-aov(genes.log2.assay.filt ~ tumor_stage, data=Genes.log2.TumorStage)
      
      ## letrinhas
      GroupLetrinhasTumorStage <- HSD.test(antumorstage, 'tumor_stage', alpha = 0.05, group=TRUE)$groups
      
      mstageI = round(mean(stageI$genes.log2.assay.filt),2)
      mstageIi = round(mean(stageIi$genes.log2.assay.filt),2)
      mstageIii = round(mean(stageIii$genes.log2.assay.filt),2)
      mstageIv = round(mean(stageIv$genes.log2.assay.filt),2)
      
      legendastageI = paste("Mean =",mstageI,"\nn =",qtdstageI)
      legendastageIi = paste("Mean =",mstageIi,"\nn =",qtdstageIi)
      legendastageIii = paste("Mean =",mstageIii,"\nn =",qtdstageIii)
      legendastageIv = paste("Mean =",mstageIv,"\nn =",qtdstageIv)
      
      
      
      p.anovaTumorStage =  as.data.frame(unlist(summary(antumorstage)))
      
      if (p.anovaTumorStage['Pr(>F)1',] < 0.001){
        legtoprightTumorStage =  paste('P = ', format(p.anovaTumorStage['Pr(>F)1',],scientific = T))
      }else{
        legtoprightTumorStage = paste("P = ",(round(p.anovaTumorStage['Pr(>F)1',],3)))
      }
      
      ### Parametric t test
      TipoTeste = "Anova with Benjamini-Hochberg correction"
    }
    
    NameOfImage = paste(NameOfGenes[i]," expression among Tumor Stage.jpg",sep="")
    
    DirectoryOfImage = paste(paste(CaminhoCompleto,NameOfGenes[i],sep="/"),NameOfImage,sep="/")
    
    png(file=DirectoryOfImage,
        width=800, height=600)  
    #tiff(file=DirectoryOfImage,res = 300)
    LabelMolecular = NameOfGenes[i]
    bx_color4 <- c("#FF6E00","brown3","lightseagreen","mediumorchid3")
    
    par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1 ## bottom, left, top, and right
    
    boxplot(genes.log2.assay.filt ~ tumor_stage, 
            data=Genes.log2.TumorStage, varwidth=F,
            #main=bquote(expression(italic(.(LabelMolecular)))), Pra baixo
            medlwd=2,
            plot=T,outline=T,outcol="white",par(font.lab=1),
            names=c("Stage I","Stage II","Stage III","Stage IV"),
            #ylim=c(9.3,13.5),
            las=1,
            boxwex=0.40,par(cex.lab=1.3),par(cex.axis=1.3),
            boxcol=bx_color4, whiskcol=bx_color4,
            medcol=bx_color4,staplecol=bx_color4,
            col=alpha(bx_color4, 0.5))
    
    axis(side=1, at=1:4, labels=c("","","",""),lty=1,las=1)
    
    title(ylab = "Relative gene expression (log2)", cex.lab = 1.3, line = 3.5)
    title(main=bquote(italic(.(NameOfGenes[i]))), cex.main=1.3, font.main=1)
    
    mtext(c(legendastageI,legendastageIi,legendastageIii,legendastageIv),
          side=1,line=3.5, at=1:4, cex=1.3)
    
    stripchart(genes.log2.assay.filt ~ tumor_stage,
               vertical=TRUE, data=Genes.log2.TumorStage, 
               method="jitter", add=TRUE, pch=19, col=alpha(bx_color5, 0.8))
    
    
    legend("topright",legend=bquote(italic(.(legtoprightTumorStage))),bty="n",cex=1.3)
    
    
    if ( table(GroupLetrinhasTumorStage[,2])[1] != 4 ){
      text(x=1, y=max(stageI$genes.log2.assay.filt), GroupLetrinhasTumorStage["CIN",2],pos=3, cex=1.3)
      text(x=2, y=max(stageIi$genes.log2.assay.filt),GroupLetrinhasTumorStage["EBV",2],pos=3, cex=1.3)
      text(x=3, y=max(stageIii$genes.log2.assay.filt),GroupLetrinhasTumorStage["GS",2],pos=3, cex=1.3)
      text(x=4, y=max(stageIv$genes.log2.assay.filt),GroupLetrinhasTumorStage["MSI",2],pos=3, cex=1.3)
    }
    
    
    dev.off()
    
    ## END Tumor Stage
  
    ## TCGA comparing SELECT GENE exp with TCGA molecular subtypes and normal samples ---------------
    
    Genes.log2.molecular <- Genes.log2.subtypes
    Genes.log2.molecular1 <- subset(Genes.log2.molecular, Genes.log2.molecular$subtype_Molecular.Subtype=="CIN")
    Genes.log2.molecular1 <- rbind(Genes.log2.molecular1, subset(Genes.log2.molecular, Genes.log2.molecular$subtype_Molecular.Subtype=="EBV"))
    Genes.log2.molecular1 <- rbind(Genes.log2.molecular1, subset(Genes.log2.molecular, Genes.log2.molecular$subtype_Molecular.Subtype=="GS"))
    Genes.log2.molecular1 <- rbind(Genes.log2.molecular1, subset(Genes.log2.molecular, Genes.log2.molecular$subtype_Molecular.Subtype=="MSI"))
    
    group_NT <- subset(Genes.log2.molecular, Genes.log2.molecular$Group=="NT")
    group_NT$subtype_Lauren.Class<-"NT"
    group_NT$subtype_Molecular.Subtype<-"NT"
    group_NT$subtype_WHO.Class<-"NT"
    
    Genes.log2.molecular1 <- rbind(Genes.log2.molecular1, subset(group_NT, group_NT$subtype_Molecular.Subtype=="NT"))
    
    #table(Genes.log2.molecular1$subtype_Molecular.Subtype)
    #CIN EBV  GS MSI  NT 
    #138  25  54  60  35  
    
    CIN = Genes.log2.molecular1[Genes.log2.molecular1$subtype_Molecular.Subtype == 'CIN',]
    EBV = Genes.log2.molecular1[Genes.log2.molecular1$subtype_Molecular.Subtype == 'EBV',]
    GS = Genes.log2.molecular1[Genes.log2.molecular1$subtype_Molecular.Subtype == 'GS',]
    MSI = Genes.log2.molecular1[Genes.log2.molecular1$subtype_Molecular.Subtype == 'MSI',]
    NTMolecular = Genes.log2.molecular1[Genes.log2.molecular1$subtype_Molecular.Subtype == 'NT',]
    qtdCIN = nrow(CIN)
    qtdEBV = nrow(EBV)
    qtdGS = nrow(GS)
    qtdMSI = nrow(MSI)
    qtdNTMolecular = nrow(NTMolecular)
    
    ## shapiro
    STCIN = shapiro.test(subset(Genes.log2.molecular1$genes.log2.assay.filt, Genes.log2.molecular1$subtype_Molecular.Subtype == "CIN"))
    STEBV = shapiro.test(subset(Genes.log2.molecular1$genes.log2.assay.filt, Genes.log2.molecular1$subtype_Molecular.Subtype == "EBV"))
    STGS  = shapiro.test(subset(Genes.log2.molecular1$genes.log2.assay.filt, Genes.log2.molecular1$subtype_Molecular.Subtype == "GS"))
    STMSI = shapiro.test(subset(Genes.log2.molecular1$genes.log2.assay.filt, Genes.log2.molecular1$subtype_Molecular.Subtype == "MSI"))
    STNT = shapiro.test(subset(Genes.log2.molecular1$genes.log2.assay.filt, Genes.log2.molecular1$subtype_Molecular.Subtype == "NT"))
    
    if ((STCIN$p.value < 0.05 ) | (STEBV$p.value < 0.05)| (STGS$p.value < 0.05)| (STMSI$p.value < 0.05)| (STNT$p.value < 0.05)){
      ## Non parametric Kruskal.Wallis Test
      kw.stad.molecular <- kruskal.test(genes.log2.assay.filt ~ subtype_Molecular.Subtype, data=Genes.log2.molecular1)
      
      mCIN = round(median(CIN$genes.log2.assay.filt),2)
      mEBV = round(median(EBV$genes.log2.assay.filt),2)
      mGS = round(median(GS$genes.log2.assay.filt),2)
      mMSI = round(median(MSI$genes.log2.assay.filt),2)
      mNTMolecular = round(median(NTMolecular$genes.log2.assay.filt),2)
      
      legendaCIN = paste("Median =",mCIN,"\nn =",qtdCIN)
      legendaEBV = paste("Median =",mEBV,"\nn =",qtdEBV)
      legendaGS = paste("Median =",mGS,"\nn =",qtdGS)
      legendaMSI = paste("Median =",mMSI,"\nn =",qtdMSI)
      legendaNTMolecular = paste("Median =",mNTMolecular,"\nn =",qtdNTMolecular)
      
  
      if (kw.stad.molecular$p.value < 0.001){
        legtoprightMolecular =  paste('P = ', format(kw.stad.molecular$p.value,scientific = T))
      }else{
        legtoprightMolecular = paste("P = ",(round(kw.stad.molecular$p.value,3)))
      }
      
      GroupLetrinhasMolecular <- kruskal(Genes.log2.molecular1$genes.log2.assay.filt,
                                         Genes.log2.molecular1$subtype_Molecular.Subtype,
                                         alpha=0.025, group=TRUE, p.adj="BH")$groups

      #Não paramétrico
      TipoTeste = "Kruskal-Wallis with Benjamini-Hochberg correction"
    }else{ 
      ### Parametric t test
      anMolecular<-aov(genes.log2.assay.filt ~ subtype_Molecular.Subtype, data=Genes.log2.molecular1)
      
      ## letrinhas
      GroupLetrinhasMolecular <- HSD.test(anMolecular, 'subtype_Molecular.Subtype', alpha = 0.05, group=TRUE)$groups
      
      
      mCIN = round(mean(CIN$genes.log2.assay.filt),2)
      mEBV = round(mean(EBV$genes.log2.assay.filt),2)
      mGS = round(mean(GS$genes.log2.assay.filt),2)
      mMSI = round(mean(MSI$genes.log2.assay.filt),2)
      mNTMolecular = mean(median(NTMolecular$genes.log2.assay.filt),2)
      
      legendaCIN = paste("Mean =",mCIN,"\nn =",qtdCIN)
      legendaEBV = paste("Mean =",mEBV,"\nn =",qtdEBV)
      legendaGS = paste("Mean =",mGS,"\nn =",qtdGS)
      legendaMSI = paste("Mean =",mMSI,"\nn =",qtdMSI)
      legendaNTMolecular = paste("Mean =",mNTMolecular,"\nn =",qtdNTMolecular)
      
      p.anovaMole =  as.data.frame(unlist(summary(anMolecular)))
      
      if (p.anovaMole['Pr(>F)1',] < 0.001){
        legtoprightMolecular =  paste('P = ', format(p.anovaMole['Pr(>F)1',],scientific = T))
      }else{
        legtoprightMolecular = paste("P = ",(round(p.anovaMole['Pr(>F)1',],3)))
      }

      ### Parametric t test
      TipoTeste = "Anova with Benjamini-Hochberg correction"
    }
    
    
    NameOfImage = paste(NameOfGenes[i]," expression among TCGA molecular-subtypes.jpg",sep="")
    
    DirectoryOfImage = paste(paste(CaminhoCompleto,NameOfGenes[i],sep="/"),NameOfImage,sep="/")
    
    png(file=DirectoryOfImage,
        width=800, height=600)  
    #tiff(file=DirectoryOfImage,res = 300)
    LabelMolecular = NameOfGenes[i]
    bx_color5 <- c("#FF6E00","brown3","lightseagreen","mediumorchid3","cornflowerblue")
    
    par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1 ## bottom, left, top, and right
    
    boxplot(genes.log2.assay.filt ~ subtype_Molecular.Subtype, 
            data=Genes.log2.molecular1, varwidth=F,
            #main=bquote(expression(italic(.(LabelMolecular)))), Pra baixo
            medlwd=2,
            plot=T,outline=T,outcol="white",par(font.lab=1),
            names=c("CIN","EBV","GS","MSI","NT"),
            #ylim=c(9.3,13.5),
            las=1,
            boxwex=0.40,par(cex.lab=1.3),par(cex.axis=1.3),
            boxcol=bx_color5, whiskcol=bx_color5,
            medcol=bx_color5,staplecol=bx_color5,
            col=alpha(bx_color5, 0.5))
    
    axis(side=1, at=1:5, labels=c("","","","",""),lty=1,las=1)
    
    title(ylab = "Relative gene expression (log2)", cex.lab = 1.3, line = 3.5)
    title(main=bquote(italic(.(NameOfGenes[i]))), cex.main=1.3, font.main=1)
    
    mtext(c(legendaCIN,legendaEBV,legendaGS,legendaMSI,legendaNTMolecular),
          side=1,line=3.5, at=1:5, cex=1.3)
    
    stripchart(genes.log2.assay.filt ~ subtype_Molecular.Subtype,
               vertical=TRUE, data=Genes.log2.molecular1, 
               method="jitter", add=TRUE, pch=19, col=alpha(bx_color5, 0.8))
    
  
    legend("topright",legend=bquote(italic(.(legtoprightMolecular))),bty="n",cex=1.3)
    
    
    if ( table(GroupLetrinhasMolecular[,2])[1] != 5 ){
      text(x=1, y=max(CIN$genes.log2.assay.filt), GroupLetrinhasMolecular["CIN",2],pos=3, cex=1.3)
      text(x=2, y=max(EBV$genes.log2.assay.filt),GroupLetrinhasMolecular["EBV",2],pos=3, cex=1.3)
      text(x=3, y=max(GS$genes.log2.assay.filt),GroupLetrinhasMolecular["GS",2],pos=3, cex=1.3)
      text(x=4, y=max(MSI$genes.log2.assay.filt),GroupLetrinhasMolecular["MSI",2],pos=3, cex=1.3)
      text(x=5, y=max(NTMolecular$genes.log2.assay.filt),GroupLetrinhasMolecular["NT",2],pos=3, cex=1.3)
    }
    
    #mtext(TipoTeste,side=1,line=4, at=2.5, cex=0.9)
    
    dev.off()
    
    ## END TCGA
    
    
    ## WHO  comparing SELECT GENE exp with WHO  and normal samples ---------------
    
    
    Genes.log2.who <- subset(Genes.log2.subtypes[1:7],subtype_WHO.Class !="NA")
    colnames(Genes.log2.who)[4:7]<-c("barcode","subtype_Lauren.Class","subtype_Molecular.Subtype","subtype_WHO.Class")
    Genes.log2.who1 <- subset(Genes.log2.who, Genes.log2.who$subtype_WHO.Class=="Mixed")
    Genes.log2.who1 <- rbind(Genes.log2.who1, subset(Genes.log2.who, Genes.log2.who$subtype_WHO.Class=="Mucinous"))
    Genes.log2.who1 <- rbind(Genes.log2.who1, subset(Genes.log2.who, Genes.log2.who$subtype_WHO.Class=="Papillary"))
    Genes.log2.who1 <- rbind(Genes.log2.who1, subset(Genes.log2.who, Genes.log2.who$subtype_WHO.Class=="Poorly_Cohesive"))
    Genes.log2.who1 <- rbind(Genes.log2.who1, subset(Genes.log2.who, Genes.log2.who$subtype_WHO.Class=="Tubular"))
    Genes.log2.who1 <- rbind(Genes.log2.who1, subset(group_NT[1:7], group_NT$subtype_WHO.Class=="NT"))
    
    
    
    
    Genes.log2.who1$subtype_WHO.Class = as.character(Genes.log2.who1$subtype_WHO.Class)
    Genes.log2.who1$subtype_WHO.Class = as.factor(Genes.log2.who1$subtype_WHO.Class)
    #table(Genes.log2.who1$subtype_WHO.Class)
    #Mixed  Mucinous  Papillary   Poorly_Cohesive  Tubular  NT
    #19     16        21          66               128      35
    
    
    WHOMixed = Genes.log2.who1[Genes.log2.who1$subtype_WHO.Class == "Mixed",]
    WHOMucinous = Genes.log2.who1[Genes.log2.who1$subtype_WHO.Class == "Mucinous",]
    WHOPapillary = Genes.log2.who1[Genes.log2.who1$subtype_WHO.Class == "Papillary",]  
    WHOPoorly_Cohesive = Genes.log2.who1[Genes.log2.who1$subtype_WHO.Class == "Poorly_Cohesive",]
    WHOTubular = Genes.log2.who1[Genes.log2.who1$subtype_WHO.Class == "Tubular",]
    WHONT =Genes.log2.who1[Genes.log2.who1$subtype_WHO.Class == "NT",]
    
    QtdWHOMixed = nrow(WHOMixed)
    QtdWHOMucinous = nrow(WHOMucinous)
    QtdWHOPapillary = nrow(WHOPapillary)
    QtdWHOPoorly_Cohesive = nrow(WHOPoorly_Cohesive)
    QtdWHOTubular = nrow(WHOTubular)
    QtdWHONT = nrow(WHONT)
    
    
    
    ## shapiro
    SWHOMixed = shapiro.test(subset(Genes.log2.who1$genes.log2.assay.filt, Genes.log2.who1$subtype_WHO.Class == "Mixed"))
    SWHOMucinous = shapiro.test(subset(Genes.log2.who1$genes.log2.assay.filt, Genes.log2.who1$subtype_WHO.Class == "Mucinous"))
    SWHOPapillary = shapiro.test(subset(Genes.log2.who1$genes.log2.assay.filt, Genes.log2.who1$subtype_WHO.Class == "Papillary"))
    SWHOPoorly_Cohesive = shapiro.test(subset(Genes.log2.who1$genes.log2.assay.filt, Genes.log2.who1$subtype_WHO.Class == "Poorly_Cohesive"))
    SWHOTubular = shapiro.test(subset(Genes.log2.who1$genes.log2.assay.filt, Genes.log2.who1$subtype_WHO.Class == "Tubular"))
    SWHONT = shapiro.test(subset(Genes.log2.who1$genes.log2.assay.filt, Genes.log2.who1$subtype_WHO.Class == "NT"))
    
    if ((SWHOMixed$p.value < 0.05 ) | (SWHOMucinous$p.value < 0.05)| (SWHOPapillary$p.value < 0.05)| (SWHOPoorly_Cohesive$p.value < 0.05)|
        (SWHOTubular$p.value < 0.05) | (SWHONT$p.value < 0.05)){
      ## Non parametric Kruskal.Wallis Test
      kw.stad.who <- kruskal.test(genes.log2.assay.filt ~ subtype_WHO.Class, data=Genes.log2.who1)

      mWHOMixed = round(median(WHOMixed$genes.log2.assay.filt),2)
      mWHOMucinous = round(median(WHOMucinous$genes.log2.assay.filt),2)
      mWHOPapillary = round(median(WHOPapillary$genes.log2.assay.filt),2)
      mWHOPoorly_Cohesive = round(median(WHOPoorly_Cohesive$genes.log2.assay.filt),2)
      mWHOTubular = round(median(WHOTubular$genes.log2.assay.filt),2)
      mWHONT = round(median(WHONT$genes.log2.assay.filt),2)
      
      WHOlegendaMixed = paste("Median =",mWHOMixed,"\nn =",QtdWHOMixed)
      WHOlegendaMucinous = paste("Median =",mWHOMucinous,"\nn =",QtdWHOMucinous)
      WHOlegendaPapillary = paste("Median =",mWHOPapillary,"\nn =",QtdWHOPapillary)
      WHOlegendaPoorly_Cohesive = paste("Median =",mWHOPoorly_Cohesive,"\nn =",QtdWHOPoorly_Cohesive)
      WHOlegendaTubular = paste("Median =",mWHOTubular,"\nn =",QtdWHOTubular)
      WHOlegendaNT = paste("Median =",mWHONT,"\nn =",QtdWHONT)
      
      
      if (kw.stad.who$p.value < 0.001){
        legtoprightWHO =  paste('P = ', format(kw.stad.who$p.value,scientific = T))
      }else{
        legtoprightWHO = paste("P = ",(round(kw.stad.who$p.value,3)))
      }
      
      GroupLetrinhasWHO <- kruskal(Genes.log2.who1$genes.log2.assay.filt,
                                         Genes.log2.who1$subtype_WHO.Class,
                                         alpha=0.025, group=TRUE, p.adj="BH")$groups
      
      #Não paramétrico
      TipoTeste = "Kruskal-Wallis with Benjamini-Hochberg correction"
      
      
    }else{ 
      ### Parametric t test
      anWHO<-aov(genes.log2.assay.filt ~ subtype_WHO.Class, data=Genes.log2.who1)
      
      ## letrinhas
      GroupLetrinhasWHO <- HSD.test(anWHO, 'subtype_WHO.Class', alpha = 0.05, group=TRUE)$groups
      
      
      mWHOMixed = round(mean(WHOMixed$genes.log2.assay.filt),2)
      mWHOMucinous = round(mean(WHOMucinous$genes.log2.assay.filt),2)
      mWHOPapillary = round(mean(WHOPapillary$genes.log2.assay.filt),2)
      mWHOPoorly_Cohesive = round(mean(WHOPoorly_Cohesive$genes.log2.assay.filt),2)
      mWHOTubular = round(mean(WHOTubular$genes.log2.assay.filt),2)
      mWHONT = round(mean(WHONT$genes.log2.assay.filt),2)
      
      
      WHOlegendaMixed = paste("Mean =",mWHOMixed,"\nn =",QtdWHOMixed)
      WHOlegendaMucinous = paste("Mean =",mWHOMucinous,"\nn =",QtdWHOMucinous)
      WHOlegendaPapillary = paste("Mean =",mWHOPapillary,"\nn =",QtdWHOPapillary)
      WHOlegendaPoorly_Cohesive = paste("Mean =",mWHOPoorly_Cohesive,"\nn =",QtdWHOPoorly_Cohesive)
      WHOlegendaTubular = paste("Mean =",mWHOTubular,"\nn =",QtdWHOTubular)
      WHOlegendaNT = paste("Mean =",mWHONT,"\nn =",QtdWHONT)

      p.anovaWHO =  as.data.frame(unlist(summary(anWHO)))
      
      if (p.anovaWHO['Pr(>F)1',] < 0.001){
        legtoprightWHO =  paste('P = ', format(p.anovaWHO['Pr(>F)1',],scientific = T))
      }else{
        legtoprightWHO = paste("P = ",(round(p.anovaWHO['Pr(>F)1',],3)))
      }
      
      ### Parametric t test
      TipoTeste = "Anova with Benjamini-Hochberg correction"
    }
    
    
    NameOfImage = paste(NameOfGenes[i]," expression among WHO classification.jpg",sep="")
    
    DirectoryOfImage = paste(paste(CaminhoCompleto,NameOfGenes[i],sep="/"),NameOfImage,sep="/")
    
    png(file=DirectoryOfImage,
        width=800, height=600)  
    #tiff(file=DirectoryOfImage,res = 300)
    LabelWHO = NameOfGenes[i]
    bx_color6 <- c("#FF6E00","lightseagreen","gold3","brown3","mediumorchid3","cornflowerblue")
    
    par(mar = c(5,5,4,0) + 0.1) ## default is c(5,4,4,2) + 0.1 ## bottom, left, top, and right
    
    
    
    boxplot(genes.log2.assay.filt ~ subtype_WHO.Class, 
            data=Genes.log2.who1, varwidth=F,
            #main=bquote(expression(italic(.(LabelMolecular)))), Pra baixo
            medlwd=2,
            plot=T,outline=T,outcol="white",par(font.lab=1),
            names=c("Mixed","Mucinous", "NT","Papillary","Poorly cohesive","Tubular"),
            #ylim=c(9.3,13.5),
            las=1,
            boxwex=0.40,par(cex.lab=1.3),par(cex.axis=1.3),
            boxcol=bx_color6, whiskcol=bx_color6,
            medcol=bx_color6,staplecol=bx_color6,
            col=alpha(bx_color6, 0.5))
  
    
    axis(side=1, at=1:6, labels=c("","","","","",""),lty=1,las=1)
    
    
    title(ylab = "Relative gene expression (log2)", cex.lab = 1.3, line = 3.5)
    title(main=bquote(italic(.(NameOfGenes[i]))), cex.main=1.3, font.main=1)
    
    mtext(c(WHOlegendaMixed,WHOlegendaMucinous,WHOlegendaNT,WHOlegendaPapillary,WHOlegendaPoorly_Cohesive,WHOlegendaTubular),
          side=1,line=3.5, at=1:6, cex=1.3)
    
    stripchart(genes.log2.assay.filt ~ subtype_WHO.Class,
               vertical=TRUE, data=Genes.log2.who1, 
               method="jitter", add=TRUE, pch=19, col=alpha(bx_color6, 0.8))
    
    
    legend("topright",legend=bquote(italic(.(legtoprightWHO))),bty="n",cex=1.3)
    
    
    
    if ( table(GroupLetrinhasWHO[,2])[1] != 6 ){
      text(x=1, y=max(WHOMixed$genes.log2.assay.filt), GroupLetrinhasWHO["Mixed",2],pos=3, cex=1.3)
      text(x=2, y=max(WHOMucinous$genes.log2.assay.filt),GroupLetrinhasWHO["Mucinous",2],pos=3, cex=1.3)
      text(x=3, y=max(WHONT$genes.log2.assay.filt),GroupLetrinhasWHO["NT",2],pos=3, cex=1.3)
      text(x=4, y=max(WHOPapillary$genes.log2.assay.filt),GroupLetrinhasWHO["Papillary",2],pos=3, cex=1.3)
      text(x=5, y=max(WHOPoorly_Cohesive$genes.log2.assay.filt),GroupLetrinhasWHO["Poorly_Cohesive",2],pos=3, cex=1.3)
      text(x=6, y=max(WHOTubular$genes.log2.assay.filt),GroupLetrinhasWHO["Tubular",2],pos=3, cex=1.3)
      
    }
    
    #mtext(TipoTeste,side=1,line=4, at=2.5, cex=0.9)
    
    dev.off()
    
    ## END WHO 
    
    
    
    ### Preparing clinical data -------
    
    stad.clinical_test <- stad.exp@colData
    stad.clinical = as.data.frame(stad.clinical_test)
    
    
    rownames(stad.assay) <- gsub( "[|].*$", "", rownames(stad.assay))
    
    stad.clinical <- cbind(stad.clinical, stad.assay[NameOfGenes[i],])
    colnames(stad.clinical)[length(stad.clinical)] <- "genes.stad.assay"
    
    # Teste para ver qual a maior data
    stad.clinical$days_to_death = as.numeric(stad.clinical$days_to_death)
    stad.clinical$days_to_last_follow_up = as.numeric(stad.clinical$days_to_last_follow_up)
    stad.clinical$overall_survival = NA
    
    stad.clinical <- stad.clinical[with(stad.clinical,order(stad.clinical$patient)),]
    
    for (t in 1:length(stad.clinical$days_to_death)) {
      
      if (is.na(stad.clinical$days_to_death[t]) & is.na(stad.clinical$days_to_last_follow_up[t]) ){
        stad.clinical$overall_survival[t] = NA
      }else if (is.na(stad.clinical$days_to_death[t])){
        stad.clinical$overall_survival[t] = stad.clinical$days_to_last_follow_up[t]
      }else if ( is.na(stad.clinical$days_to_last_follow_up[t])){
        stad.clinical$overall_survival[t] = stad.clinical$days_to_death[t]
      }else if ( stad.clinical$days_to_death[t] > stad.clinical$days_to_last_follow_up[t]){
        stad.clinical$overall_survival[t] = stad.clinical$days_to_death[t]
      }else if ( stad.clinical$days_to_death[t] <= stad.clinical$days_to_last_follow_up[t]){
        stad.clinical$overall_survival[t] = stad.clinical$days_to_last_follow_up[t]
      }
      
    }
    
    stad.clinical$genes.stad.assay <- as.numeric(stad.clinical$genes.stad.assay)
    ## transform days into months = STADsurvival_info$m_overall_survival
    stad.clinical$m_overall_survival <- stad.clinical$overall_survival/30
    stad.clinical$m_overall_survival <- as.numeric(stad.clinical$m_overall_survival)
    # recode vital status to 0=alive 1=dead
    stad.clinical$vital_status1 <- stad.clinical$vital_status
    stad.clinical$vital_status1[stad.clinical$vital_status=="alive"] <- 0
    stad.clinical$vital_status1[stad.clinical$vital_status=="dead"] <- 1
    stad.clinical$vital_status1 <- as.numeric(stad.clinical$vital_status1)
    
    
    ## testing if there is some difference between TULP3 and stomach cancer subtypes
    #which(colnames(stad.clinical)=="subtype_Molecular.Subtype")
    km_genes_lauren <- stad.clinical[,c(grep("barcode", colnames(stad.clinical)),grep("m_overall_survival", colnames(stad.clinical)),
                                        grep("vital_status1", colnames(stad.clinical)),
                                        grep("subtype_Lauren.Class", colnames(stad.clinical)),grep("genes.stad.assay", colnames(stad.clinical)))]
    
    km_genes_lauren$subtype_Lauren.Class = as.character(km_genes_lauren$subtype_Lauren.Class)
    km_genes_lauren = subset(km_genes_lauren,km_genes_lauren$subtype_Lauren.Class != "NA")
    km_genes_lauren$subtype_Lauren.Class = as.factor(km_genes_lauren$subtype_Lauren.Class)
    qtdGroup =   table(km_genes_lauren$subtype_Lauren.Class)
    #qtdGroup
    #Diffuse  Intestinal  Mixed 
    #66       181         19 
    
    # Diffuse Survival probability.jpg --------------------------------
    
    ## Diffuse type
    km_genes_diffuse <- subset(km_genes_lauren, km_genes_lauren$subtype_Lauren.Class == "Diffuse")
    km_genes_diffuse <- km_genes_diffuse[with(km_genes_diffuse,order(km_genes_diffuse$genes.stad.assay)),]
    ## remove NA values
    km_genes_diffuse<-na.omit(km_genes_diffuse)
    
    rownames(km_genes_diffuse)<-c(1: as.numeric(qtdGroup["Diffuse"]))
    
    lows = subset(km_genes_diffuse,km_genes_diffuse$genes.stad.assay <= median(km_genes_diffuse$genes.stad.assay))
    high = subset(km_genes_diffuse,km_genes_diffuse$genes.stad.assay > median(km_genes_diffuse$genes.stad.assay))
    hRatio = hazard.ratio(c(rep(0,as.numeric(lengths(lows)[1])),
                            rep(1,as.numeric(lengths(high)[1]))),
                          km_genes_diffuse$m_overall_survival,km_genes_diffuse$vital_status1)
    diffuse <- survfit(Surv(km_genes_diffuse$m_overall_survival,km_genes_diffuse$vital_status1) ~ 
                         c(rep(0,as.numeric(lengths(lows)[1])),
                           rep(1,as.numeric(lengths(high)[1]))),
                       conf.int=0.95)
    
    lowMedian = round(as.numeric(summary(diffuse)$table[, "median"][1]),2) #low
    highMedian = round(as.numeric(summary(diffuse)$table[, "median"][2]),2) #high
    
    #Tratamento para a quantidade de meses
    maiorSobrevida = max(km_genes_diffuse$m_overall_survival)
    maiorSobrevida = (1 + trunc(maiorSobrevida/6)) * 6
    resultDiv =  maiorSobrevida / 6
    
    times  = c(0)
    diffuse_low = c(as.numeric(table(lows$m_overall_survival >= 0)["TRUE"]))
    diffuse_high = diffuse_high = c(as.numeric(table(high$m_overall_survival >= 0)["TRUE"]))
    for (j in  1:resultDiv){
      times = c(times,j*6)
      if (is.na(as.numeric(table(lows$m_overall_survival >= j*6)["TRUE"]))){
        diffuse_low = c(diffuse_low,0)
      }else{
        diffuse_low = c(diffuse_low,as.numeric(table(lows$m_overall_survival >= j*6)["TRUE"])) 
      }
      
      if (is.na(as.numeric(table(high$m_overall_survival >= j*6)["TRUE"]))){
        diffuse_high = c(diffuse_high,0)
      }else{
        diffuse_high = c(diffuse_high,as.numeric(table(high$m_overall_survival >= j*6)["TRUE"]))
      }  
      
    }
    
    legendTopRight = c(paste("HR = ",round(hRatio$hazard.ratio,2)," (",
                             round(hRatio$lower,2),"-",round(hRatio$upper,2),")",sep=""),paste("P = ",round(hRatio$p.value,3),sep=""))
    
    
    ##Colocar aqui a mediana do low e do high  ###SEMPRE A MEDIANA
    legendofcut = c(paste("Median survival group low =", lowMedian),
                    paste("Median survival group high = ", highMedian))
    
    # number at risk
    summaryTeste = summary(diffuse, times)
    
    NameOfImage = paste(NameOfGenes[i]," Survival probability of Diffuse-type.jpg",sep="")
    
    DirectoryOfImage = paste(paste(CaminhoCompleto,NameOfGenes[i],sep="/"),NameOfImage,sep="/")
    
    png(file=DirectoryOfImage,
        width=800, height=600)  
    #tiff(file=DirectoryOfImage,res = 300)
    par(mar = c(9,5,4,1.5) + 0.1) ## (9,5,4,1.5) default is c(5,4,4,2) + 0.1 ## bottom, left, top, and right
    
    plot(diffuse, col=c('blue','red2'), lty=1, lwd=2, cex.lab=1.3, cex.axis=1.3, las=1, xaxt="n", mark.time=T)
    
    title(ylab = "Diffuse-type survival probability", cex.lab = 1.3, line = 3.0)
    
    axis(side=1, at=times, labels=c(times), cex.axis=1.2, las=0.8)
    
    #abline(h=0.5, col="gray40", lty=2)
    title(main=bquote(italic(.(NameOfGenes[i]))),cex.main=1.3, font.main=1)
    legend("bottomleft",
           legend=legendofcut,
           bty="n", cex=1.3, col=c('blue','red2'), lty=1, lwd=2.5, y.intersp=1)
    legend("topright", cex=1.3, bty="n", y.intersp=1.3,
           legend=legendTopRight)
    
    mtext("Number at risk", side=1, line=2.8, cex=1.3, at=-1.15)
    mtext("Low", side=1, line=4, cex=1.3, at=-6.8, col="blue")
    mtext(diffuse_low, side=1, line=4, cex=1.3, at=c(times))
    mtext("High", side=1, line=5.3, cex=1.3, at=-6.6, col="red2")
    mtext(diffuse_high, side=1, line=5.3, cex=1.3, at=c(times))
    ##
    dev.off()
    
    ## END Survival probability of Diffuse-type.jpg    
    
    # Intestinal Survival probability ---------------------------------
    ## Intestinal type
    km_genes_intestinal <- subset(km_genes_lauren, km_genes_lauren$subtype_Lauren.Class == "Intestinal")
    km_genes_intestinal<-na.omit(km_genes_intestinal)
    km_genes_intestinal <- km_genes_intestinal[with(km_genes_intestinal,order(km_genes_intestinal$genes.stad.assay)),]
    rownames(km_genes_intestinal)<-c(1: as.numeric(lengths(km_genes_intestinal)[1]))
    
    intlows = subset(km_genes_intestinal,km_genes_intestinal$genes.stad.assay <= median(km_genes_intestinal$genes.stad.assay))
    inthigh = subset(km_genes_intestinal,km_genes_intestinal$genes.stad.assay > median(km_genes_intestinal$genes.stad.assay))
    
    
    hratioInt =  hazard.ratio(c(rep(0,as.numeric(lengths(intlows)[1])),
                                rep(1,as.numeric(lengths(inthigh)[1]))),
                              km_genes_intestinal$m_overall_survival,km_genes_intestinal$vital_status1)
    
    
    intestinal <- survfit(Surv(km_genes_intestinal$m_overall_survival,km_genes_intestinal$vital_status1) ~ 
                            c(rep(0,as.numeric(lengths(intlows)[1])),
                              rep(1,as.numeric(lengths(inthigh)[1]))),
                          conf.int=0.95)
    
    
    lowMedianInt = round(as.numeric(summary(intestinal)$table[, "median"][1]),2) #low
    highMedianInt = round(as.numeric(summary(intestinal)$table[, "median"][2]),2) #high
    
    #Tratamento para a quantidade de meses
    maiorSobrevidaInt = max(km_genes_intestinal$m_overall_survival)
    maiorSobrevidaInt = (1 + trunc(maiorSobrevidaInt/6)) * 6
    resultDivInt =  maiorSobrevidaInt / 6
    
    timesIntestinal  = c(0)
    intestinal_low = c(as.numeric(table(intlows$m_overall_survival >= 0)["TRUE"]))
    intestinal_high  = c(as.numeric(table(inthigh$m_overall_survival >= 0)["TRUE"]))
    for (k in  1:resultDivInt){
      timesIntestinal = c(timesIntestinal,k*6)
      if (is.na(as.numeric(table(intlows$m_overall_survival >= k*6)["TRUE"]))){
        intestinal_low = c(intestinal_low,0)
      }else{
        intestinal_low = c(intestinal_low,as.numeric(table(intlows$m_overall_survival >= k*6)["TRUE"])) 
      }
      
      if (is.na(as.numeric(table(inthigh$m_overall_survival >= k*6)["TRUE"]))){
        intestinal_high = c(intestinal_high,0)
      }else{
        intestinal_high = c(intestinal_high,as.numeric(table(inthigh$m_overall_survival >= k*6)["TRUE"]))
      }  
      
    }
    
    legendTopRightIntestinal = c(paste("HR = ",round(hratioInt$hazard.ratio,2)," (",
                                       round(hratioInt$lower,2),"-",round(hratioInt$upper,2),")",sep=""),paste("P = ",round(hratioInt$p.value,3),sep=""))
    
    
    ###Colocar a mediana #
    legendofcutIntestinal = c(paste("Median survival group low =", lowMedianInt),
                              paste("Median survival group high = ", highMedianInt))
    
    #summary(intestinal, timesIntestinal)
    NameOfImage = paste(NameOfGenes[i]," Survival probability of Intestinal-type.jpg",sep="")
    
    DirectoryOfImage = paste(paste(CaminhoCompleto,NameOfGenes[i],sep="/"),NameOfImage,sep="/")
    
    png(file=DirectoryOfImage,
        width=800, height=600)  
    #tiff(file=DirectoryOfImage,res = 300)
    par(mar = c(9,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1 ## bottom, left, top, and right
    
    plot(intestinal, col=c('blue','red2'), lty=1, lwd=2, cex.lab=1.3, cex.axis=1.3, las=1, xaxt="n", mark.time=T)
    title(ylab = "Intestinal-type survival probability", cex.lab = 1.3, line = 3.5)
    
    axis(side=1, at=timesIntestinal,labels=timesIntestinal, cex.axis=1.2, las=1)
    #abline(h=0.5, col="gray40", lty=2)
    title(main=bquote(italic(.(NameOfGenes[i]))), cex.main=1.3, font.main=1)
    
    mtext("Number at risk", side=1, line=2.8, cex=1.3, at=-1.9)
    mtext("Low", side=1, line=4, cex=1.3, at=-6, col="blue")
    mtext(intestinal_low, side=1, line=4, cex=1.3, at=timesIntestinal)
    mtext("High", side=1, line=5.3, cex=1.3, at=-5.9, col="red2")
    mtext(intestinal_high, side=1, line=5.3, cex=1.3, at=timesIntestinal)
    
    legend("bottomleft",
           legend=legendofcutIntestinal,
           bty="n", cex=1.3, col=c('blue','red2'), lty=1, lwd=2.5, y.intersp=1)
    legend("topright", cex=1.3, bty="n", y.intersp=1.3,
           legend=legendTopRightIntestinal)
    ##
    dev.off()
    
    ## END Survival probability of Intestinal-type    
    
    
    
    # Int X Dif Survival probability between Intestinal-type and Diffuse-type -----------
    km_genes_dif_int <- subset(km_genes_lauren, km_genes_lauren$subtype_Lauren.Class == "Intestinal")
    km_genes_dif_int <- rbind(km_genes_dif_int, subset(km_genes_lauren, km_genes_lauren$subtype_Lauren.Class == "Diffuse"))
    km_genes_dif_int<-na.omit(km_genes_dif_int)
    
    qtd_Dif = subset(km_genes_dif_int,km_genes_dif_int$subtype_Lauren.Class == "Diffuse")
    qtd_Int = subset(km_genes_dif_int,km_genes_dif_int$subtype_Lauren.Class == "Intestinal")
    
    hratioDifInt =  hazard.ratio(c(rep(0,as.numeric(lengths(qtd_Int)[1])),
                                   rep(1,as.numeric(lengths(qtd_Dif)[1]))),
                                 km_genes_dif_int$m_overall_survival,km_genes_dif_int$vital_status1)
    
    dif_int <- survfit(Surv(km_genes_dif_int$m_overall_survival,km_genes_dif_int$vital_status1) ~ 
                         c(rep(0,as.numeric(lengths(qtd_Int)[1])),
                           rep(1,as.numeric(lengths(qtd_Dif)[1]))),
                       conf.int=0.95)
    
    MedianInt = round(as.numeric(summary(dif_int)$table[, "median"][1]),2) #low
    MedianDif = round(as.numeric(summary(dif_int)$table[, "median"][2]),2) #high
    
    #Tratamento para a quantidade de meses
    maiorSobrevidaDifInt = max(km_genes_dif_int$m_overall_survival)
    maiorSobrevidaDifInt = (1 + trunc(maiorSobrevidaDifInt/6)) * 6
    resultDivDifInt =  maiorSobrevidaDifInt / 6
    
    timesDifInt  = c(0)
    #Low
    dif_int_INT = c(as.numeric(table(qtd_Int$m_overall_survival >= 0)["TRUE"]))
    #High
    dif_int_DIF  = c(as.numeric(table(qtd_Dif$m_overall_survival >= 0)["TRUE"]))
    for (x in  1:resultDivDifInt){
      timesDifInt = c(timesDifInt,x*6)
      if (is.na(as.numeric(table(qtd_Int$m_overall_survival >= x*6)["TRUE"]))){
        dif_int_INT = c(dif_int_INT,0)
      }else{
        dif_int_INT = c(dif_int_INT,as.numeric(table(qtd_Int$m_overall_survival >= x*6)["TRUE"])) 
      }
      
      if (is.na(as.numeric(table(qtd_Dif$m_overall_survival >= x*6)["TRUE"]))){
        dif_int_DIF = c(dif_int_DIF,0)
      }else{
        dif_int_DIF = c(dif_int_DIF,as.numeric(table(qtd_Dif$m_overall_survival >= x*6)["TRUE"]))
      }  
      
    }
    
    legendTopRightDifInt = c(paste("HR = ",round(hratioDifInt$hazard.ratio,2)," (",
                                   round(hratioDifInt$lower,2),"-",round(hratioDifInt$upper,2),")",sep=""),paste("P = ",round(hratioDifInt$p.value,3),sep=""))
    
    
    legendofcutDifInt = c(paste("Median survival group Intestinal =", MedianInt),
                          paste("Median survival group Diffuse = ", MedianDif))
    
    # number at risk
    NameOfImage = paste(NameOfGenes[i]," Survival probability between Intestinal-type and Diffuse-type.jpg",sep="")
    
    DirectoryOfImage = paste(paste(CaminhoCompleto,NameOfGenes[i],sep="/"),NameOfImage,sep="/")
    
    png(file=DirectoryOfImage,
        width=800, height=600)
    
    #tiff(file=DirectoryOfImage,res = 300)
    
    par(mar = c(9,6,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1 ## bottom, left, top, and right
    
    
    plot(dif_int, col=c('blue','red2'), lty=1, lwd=2, cex.lab=1.3, cex.axis=1.3, las=1, xaxt="n", mark.time=T)
    title(ylab = "Survival probability", cex.lab = 1.3, line = 4.2)
    title(main=bquote(italic(.("Intestinal-type and Diffuse-type"))), cex.main=1.3, font.main=1)
    
    axis(side=1, at=timesDifInt,labels=timesDifInt,cex.axis=1.2, las=1)
    
    mtext("Number at risk", side=1, line=2.8, cex=1.3, at=-4.2)
    mtext("Intestinal", side=1, line=4, cex=1.3, at=-7.9, col="blue")
    mtext(dif_int_INT, side=1, line=4, cex=1.3, at=timesDifInt)
    mtext("Diffuse", side=1, line=5.3, cex=1.3, at=-8.9, col="red2")
    mtext(dif_int_DIF, side=1, line=5.3, cex=1.3, at=timesDifInt)
    
    legend("topright", cex=1.3, bty="n", y.intersp=1.3,
           legend=legendTopRightDifInt)
    
    legend("bottomleft",
           legend=legendofcutDifInt,
           bty="n", cex=1.3, col=c('blue','red2'), lty=1, lwd=2.5, y.intersp=1)
    
    
    ##
    dev.off()
    
    
    
    ## END Survival probability between Intestinal-type and Diffuse-type
    
    
    # TCGA survival according STAD TCGA molecular subtypes -------------------------
    
    km_genes_molecular <- stad.clinical[,c(grep("barcode", colnames(stad.clinical)),grep("m_overall_survival", colnames(stad.clinical)),
                                           grep("vital_status1", colnames(stad.clinical)),
                                           grep("subtype_Molecular.Subtype", colnames(stad.clinical)),grep("genes.stad.assay", colnames(stad.clinical)))]
    
    ## molecular subtypes
    km_genes_molecular2 <- subset(km_genes_molecular, km_genes_molecular$subtype_Molecular.Subtype == "CIN")
    km_genes_molecular2 <- rbind(km_genes_molecular2, subset(km_genes_molecular, km_genes_molecular$subtype_Molecular.Subtype == "EBV"))
    km_genes_molecular2 <- rbind(km_genes_molecular2, subset(km_genes_molecular, km_genes_molecular$subtype_Molecular.Subtype == "GS"))
    km_genes_molecular2 <- rbind(km_genes_molecular2, subset(km_genes_molecular, km_genes_molecular$subtype_Molecular.Subtype == "MSI"))
    
    rownames(km_genes_molecular2)<-c(1:as.numeric(lengths(km_genes_molecular2)[1]))
    
    ## Remove NA values
    km_genes_molecular2 <- na.omit(km_genes_molecular2)
    
    qtdKmSubMolecular = table(km_genes_molecular2$subtype_Molecular.Subtype)
    #CIN EBV GS MSI 
    #138 25  54  60
    #136 no meu
    
    ## Molecular subtypes
    MolecularCin = subset(km_genes_molecular2,km_genes_molecular2$subtype_Molecular.Subtype == "CIN")
    MolecularEBV = subset(km_genes_molecular2,km_genes_molecular2$subtype_Molecular.Subtype == "EBV")
    MolecularGS = subset(km_genes_molecular2,km_genes_molecular2$subtype_Molecular.Subtype == "GS")
    MolecularMSI = subset(km_genes_molecular2,km_genes_molecular2$subtype_Molecular.Subtype == "MSI")
    
    ## test HR
    hrmolecular = hazard.ratio(c(rep(0,as.numeric(qtdKmSubMolecular["CIN"])),
                                 rep(1,as.numeric(qtdKmSubMolecular["EBV"])),
                                 rep(2,as.numeric(qtdKmSubMolecular["GS"])),
                                 rep(3,as.numeric(qtdKmSubMolecular["MSI"]))),
                               km_genes_molecular2$m_overall_survival,km_genes_molecular2$vital_status1)
    
    molecular <- survfit(Surv(km_genes_molecular2$m_overall_survival,km_genes_molecular2$vital_status1) ~ 
                           c(rep(0,as.numeric(qtdKmSubMolecular["CIN"])),
                             rep(1,as.numeric(qtdKmSubMolecular["EBV"])),
                             rep(2,as.numeric(qtdKmSubMolecular["GS"])),
                             rep(3,as.numeric(qtdKmSubMolecular["MSI"]))),
                         conf.int=0.95)
    
    maiorSobrevidaMolecular = max(km_genes_molecular2$m_overall_survival)
    maiorSobrevidaMolecular = (1 + trunc(maiorSobrevidaMolecular/6)) * 6
    resultDivMolecular =  maiorSobrevidaMolecular / 6
    
    timesMolecular  = c(0)
    molecular_CIN = c(as.numeric(table(MolecularCin$m_overall_survival >= 0)["TRUE"]))
    molecular_EBV = c(as.numeric(table(MolecularEBV$m_overall_survival >= 0)["TRUE"]))
    molecular_GS = c(as.numeric(table(MolecularGS$m_overall_survival >= 0)["TRUE"]))
    molecular_MSI = c(as.numeric(table(MolecularMSI$m_overall_survival >= 0)["TRUE"]))
    for (s in  1:resultDivMolecular){
      timesMolecular = c(timesMolecular,s*6)
      
      if (is.na(as.numeric(table(MolecularCin$m_overall_survival >= s*6)["TRUE"]))){
        molecular_CIN = c(molecular_CIN,0)
      }else{
        molecular_CIN = c(molecular_CIN,as.numeric(table(MolecularCin$m_overall_survival >= s*6)["TRUE"])) 
      }
      
      if (is.na(as.numeric(table(MolecularEBV$m_overall_survival >= s*6)["TRUE"]))){
        molecular_EBV = c(molecular_EBV,0)
      }else{
        molecular_EBV = c(molecular_EBV,as.numeric(table(MolecularEBV$m_overall_survival >= s*6)["TRUE"]))
      }  
      
      if (is.na(as.numeric(table(MolecularGS$m_overall_survival >= s*6)["TRUE"]))){
        molecular_GS = c(molecular_GS,0)
      }else{
        molecular_GS = c(molecular_GS,as.numeric(table(MolecularGS$m_overall_survival >= s*6)["TRUE"]))
      }
      
      if (is.na(as.numeric(table(MolecularMSI$m_overall_survival >= s*6)["TRUE"]))){
        molecular_MSI = c(molecular_MSI,0)
      }else{
        molecular_MSI = c(molecular_MSI,as.numeric(table(MolecularMSI$m_overall_survival >= s*6)["TRUE"]))
      }
    }
    
    legendTopRightMolecular = c(paste("HR = ",round(hrmolecular$hazard.ratio,2)," (",
                                      round(hrmolecular$lower,2),"-",round(hrmolecular$upper,2),")",sep=""),paste("P = ",round(hrmolecular$p.value,3),sep=""))
    
    # number at risk
    obj = summary(molecular, timesMolecular)
    
    ## Preciso ajustar os titulos conforme paper... 
    NameOfImage = paste(NameOfGenes[i]," Survival probability TCGA molecular-subtypes.jpg",sep="")
    
    DirectoryOfImage = paste(paste(CaminhoCompleto,NameOfGenes[i],sep="/"),NameOfImage,sep="/")
    
    png(file=DirectoryOfImage,
        width=800, height=600)
    #tiff(file=DirectoryOfImage,res = 300)
    par(mar = c(9,6,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1 ## bottom, left, top, and right
    
    plot(molecular, col=c('#FF6E00','#A50026','darkolivegreen4','#74ADD1'), lty=1, lwd=2, cex.lab=1.3, cex.axis=1.3, las=1,
         xaxt="n", mark.time=T)
    title(ylab = "Survival probability ", cex.lab = 1.3, line = 3.5)
    axis(side=1, at=timesMolecular,labels=timesMolecular,cex.axis=1.2, las=1)
    
    mtext("Number at risk", side=1, line=2.8, cex=1.3, at=-0.3)
    mtext("CIN", side=1, line=4, cex=1.3, at=-6.6, col="#FF6E00")
    mtext(molecular_CIN, side=1, line=4, cex=1.3, at=timesMolecular)
    mtext("EBV", side=1, line=5.3, cex=1.3, at=-5.95, col="#A50026")
    mtext(molecular_EBV, side=1, line=5.3, cex=1.3, at=timesMolecular)
    mtext("GS", side=1, line=6.6, cex=1.3, at=-6.9, col="darkolivegreen4")
    mtext(molecular_GS, side=1, line=6.6, cex=1.3, at=timesMolecular)
    mtext("MSI", side=1, line=7.9, cex=1.3, at=-6.4, col="#74ADD1")
    mtext(molecular_MSI, side=1, line=7.9, cex=1.3, at=timesMolecular)
    
    title(main=bquote(italic(.("TCGA molecular-subtypes"))), cex.main=1.3, font.main=1)
    
    legend("topright", cex=1.3, bty="n", y.intersp=1.3,
           legend=legendTopRightMolecular)
    
    dev.off()
    
    ## END survival according STAD TCGA molecular subtypes
    
    
    # WHO survival according STAD TCGA WHO classification -------------------------
    
    km_tulp3_who <- stad.clinical[,c(grep("barcode", colnames(stad.clinical)),grep("m_overall_survival", colnames(stad.clinical)),
                                     grep("vital_status1", colnames(stad.clinical)),
                                     grep("subtype_WHO.Class", colnames(stad.clinical)),grep("genes.stad.assay", colnames(stad.clinical)))]
    
    km_tulp3_who2 <- subset(km_tulp3_who, km_tulp3_who$subtype_WHO.Class == "Mixed")
    km_tulp3_who2 <- rbind(km_tulp3_who2, subset(km_tulp3_who, km_tulp3_who$subtype_WHO.Class == "Mucinous"))
    km_tulp3_who2 <- rbind(km_tulp3_who2, subset(km_tulp3_who, km_tulp3_who$subtype_WHO.Class == "Papillary"))
    km_tulp3_who2 <- rbind(km_tulp3_who2, subset(km_tulp3_who, km_tulp3_who$subtype_WHO.Class == "Poorly_Cohesive"))
    km_tulp3_who2 <- rbind(km_tulp3_who2, subset(km_tulp3_who, km_tulp3_who$subtype_WHO.Class == "Tubular"))
    
    rownames(km_tulp3_who2)<-c(1:as.numeric(lengths(km_tulp3_who2)[1]))
    
    ## Remove NA values
    km_tulp3_who2 <- na.omit(km_tulp3_who2)
    km_tulp3_who2$subtype_WHO.Class = as.character(km_tulp3_who2$subtype_WHO.Class)
    km_tulp3_who2$subtype_WHO.Class = as.factor(km_tulp3_who2$subtype_WHO.Class)
    
    WHOMixed = subset(km_tulp3_who2,km_tulp3_who2$subtype_WHO.Class == "Mixed")
    WHOMucinous = subset(km_tulp3_who2,km_tulp3_who2$subtype_WHO.Class == "Mucinous")
    WHOPapillary = subset(km_tulp3_who2,km_tulp3_who2$subtype_WHO.Class == "Papillary")
    WHOPoorly_Cohesive = subset(km_tulp3_who2,km_tulp3_who2$subtype_WHO.Class == "Poorly_Cohesive")
    WHOTubular = subset(km_tulp3_who2,km_tulp3_who2$subtype_WHO.Class == "Tubular")
    
    hratioWHO = hazard.ratio(c(rep(0,as.numeric(table(km_tulp3_who2$subtype_WHO.Class)["Mixed"])),
                               rep(1,as.numeric(table(km_tulp3_who2$subtype_WHO.Class)["Mucinous"])),
                               rep(2,as.numeric(table(km_tulp3_who2$subtype_WHO.Class)["Papillary"])),
                               rep(3,as.numeric(table(km_tulp3_who2$subtype_WHO.Class)["Poorly_Cohesive"])),
                               rep(4,as.numeric(table(km_tulp3_who2$subtype_WHO.Class)["Tubular"]))),
                             km_tulp3_who2$m_overall_survival,km_tulp3_who2$vital_status1)
    
    who <- survfit(Surv(km_tulp3_who2$m_overall_survival,km_tulp3_who2$vital_status1) ~ 
                     c(rep(0,as.numeric(table(km_tulp3_who2$subtype_WHO.Class)["Mixed"])),
                       rep(1,as.numeric(table(km_tulp3_who2$subtype_WHO.Class)["Mucinous"])),
                       rep(2,as.numeric(table(km_tulp3_who2$subtype_WHO.Class)["Papillary"])),
                       rep(3,as.numeric(table(km_tulp3_who2$subtype_WHO.Class)["Poorly_Cohesive"])),
                       rep(4,as.numeric(table(km_tulp3_who2$subtype_WHO.Class)["Tubular"]))),
                   conf.int=0.95)
    
    maiorSobrevidaWHO = max(km_genes_molecular2$m_overall_survival)
    maiorSobrevidaWHO = (1 + trunc(maiorSobrevidaWHO/6)) * 6
    resultDivWHO =  maiorSobrevidaWHO / 6
    
    timesWHO  = c(0)
    who_mixed = c(as.numeric(table(WHOMixed$m_overall_survival >= 0)["TRUE"]))
    who_mucinous = c(as.numeric(table(WHOMucinous$m_overall_survival >= 0)["TRUE"]))
    who_papillary = c(as.numeric(table(WHOPapillary$m_overall_survival >= 0)["TRUE"]))
    who_poorly = c(as.numeric(table(WHOPoorly_Cohesive$m_overall_survival >= 0)["TRUE"]))
    who_tubular = c(as.numeric(table(WHOTubular$m_overall_survival >= 0)["TRUE"]))
    
    for (z in  1:resultDivWHO){
      timesWHO = c(timesWHO,z*6)
      
      if (is.na(as.numeric(table(WHOMixed$m_overall_survival >= z*6)["TRUE"]))){
        who_mixed = c(who_mixed,0)
      }else{
        who_mixed = c(who_mixed,as.numeric(table(WHOMixed$m_overall_survival >= z*6)["TRUE"])) 
      }
      
      if (is.na(as.numeric(table(WHOMucinous$m_overall_survival >= z*6)["TRUE"]))){
        who_mucinous = c(who_mucinous,0)
      }else{
        who_mucinous = c(who_mucinous,as.numeric(table(WHOMucinous$m_overall_survival >= z*6)["TRUE"]))
      }  
      
      if (is.na(as.numeric(table(WHOPapillary$m_overall_survival >= z*6)["TRUE"]))){
        who_papillary = c(who_papillary,0)
      }else{
        who_papillary = c(who_papillary,as.numeric(table(WHOPapillary$m_overall_survival >= z*6)["TRUE"]))
      }
      
      if (is.na(as.numeric(table(WHOPoorly_Cohesive$m_overall_survival >= z*6)["TRUE"]))){
        who_poorly = c(who_poorly,0)
      }else{
        who_poorly = c(who_poorly,as.numeric(table(WHOPoorly_Cohesive$m_overall_survival >= z*6)["TRUE"]))
      }
      
      if (is.na(as.numeric(table(WHOTubular$m_overall_survival >= z*6)["TRUE"]))){
        who_tubular = c(who_tubular,0)
      }else{
        who_tubular = c(who_tubular,as.numeric(table(WHOTubular$m_overall_survival >= z*6)["TRUE"]))
      }
    }
    
    legendTopRightWHO = c(paste("HR = ",round(hratioWHO$hazard.ratio,2)," (",
                                round(hratioWHO$lower,2),"-",round(hratioWHO$upper),")",sep=""),paste("P = ",round(hratioWHO$p.value,3),sep=""))
    
    
    NameOfImage = paste(NameOfGenes[i]," Survival probability WHO classification.jpg",sep="")
    
    DirectoryOfImage = paste(paste(CaminhoCompleto,NameOfGenes[i],sep="/"),NameOfImage,sep="/")
    
    png(file=DirectoryOfImage,
        width=750, height=600)
    #tiff(file=DirectoryOfImage,res = 300)
    par(mar = c(10.3,5.5,4,1.5) + 0.1) ## default is c(5,4,4,2) + 0.1 ## bottom, left, top, and right
    
    plot(who, col=c('#FF6E00','#A50026','darkolivegreen4','#74ADD1','mediumpurple3'), lty=1, lwd=2, cex.lab=1.3, cex.axis=1.3, las=1,
         xaxt="n", mark.time=T)
    
    title(ylab = "Survival probability", cex.lab = 1.3, line = 3.5)
    
    axis(side=1, at=timesWHO,labels=timesWHO,cex.axis=1.2, las=1)
    
    mtext("Number at risk", side=1, line=2.8, cex=1.3, at=-4.5)
    mtext("Mixed", side=1, line=4, cex=1.3, at=-9.6, col="#FF6E00")
    mtext(who_mixed, side=1, line=4, cex=1.3, at=timesWHO)
    mtext("Mucinous", side=1, line=5.3, cex=1.3, at=-7.6, col="#A50026")
    mtext(who_mucinous, side=1, line=5.3, cex=1.3, at=timesWHO)
    mtext("Papillary", side=1, line=6.6, cex=1.3, at=-8, col="darkolivegreen4")
    mtext(who_papillary, side=1, line=6.6, cex=1.3, at=timesWHO)
    mtext("Poorly", side=1, line=7.9, cex=1.3, at=-9.5, col="#74ADD1")
    mtext(who_poorly, side=1, line=7.9, cex=1.3, at=timesWHO)
    mtext("Tubular", side=1, line=9.1, cex=1.3, at=-8.6, col="mediumpurple4")
    mtext(who_tubular, side=1, line=9.1, cex=1.3, at=timesWHO)
    
    title(main=bquote(italic(.("WHO subtypes"))), cex.main=1.3, font.main=1)
    legend("topright", cex=1.3, bty="n", y.intersp=1.3, legend=legendTopRightWHO)
    #
    dev.off()
    
    ## END survival according STAD TCGA WHO classification
  }else{
    print(paste("Gene",NameOfGenes[i],"not exists in data"))
  }
  
}
print("Analysis completed")