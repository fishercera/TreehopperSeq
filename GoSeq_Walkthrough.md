Working with GOSeq to do GO Term Analysis
================

This is a step-by-step walkthrough of how to perform GO term enrichment
analysis on a non-model organism using the R package goseq from
Bioconductor.

Important pre goseq steps – you need a gene-to-GOterm map which is a
list of gene identifiers and an associated GO identifier. If you produce
it with my python script, GoTermsMap.py, it looks like
this:

``` 
 ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0006807-nitrogen compound metabolic process(L=2)
 ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0007154-cell communication(L=2)
 ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0009058-biosynthetic process(L=2)
 ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0044237-cellular metabolic process(L=2)
 ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0044238-primary metabolic process(L=2)
 ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0044700-single organism signaling(L=2)
 ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0044710-single-organism metabolic process(L=2)
 ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0044763-single-organism cellular process(L=2)
 ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0050789-regulation of biological process(L=2)
 ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0051716-cellular response to stimulus(L=2)
 ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0071704-organic substance metabolic process(L=2)
 ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0043227-membrane-bounded organelle(L=2)
 ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0044421-extracellular region part(L=2)
```

but will need to be edited so that it looks more like
    this:

    ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0006807-nitrogen compound metabolic process(L=2) GO:0006807
    ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0007154-cell communication(L=2)                  GO:0007154
    ECEF_Abd_TRINITY_DN10037_c0_g1  GO:0009058-biosynthetic process(L=2)                GO:0009058

That is to say, goseq wants just the GO identifier and not the verbose
category. If you are feeling clever, you can just fix it in R with
regular expressions.

You will also need a gene id to gene length map. I decided the best way
to get this was to go to one of my sample’s genes.results output
produced from the Trinity RSEM/edgeR DE Analysis scripts. I did `cut
-f 1,4` to get the gene identifier and the effective length out of a
genes.results file.

What you input to R should look like this:

    gene_id                         length
    ECEF_Abd_TRINITY_DN10007_c0_g1  1148.00
    ECEF_Abd_TRINITY_DN10023_c0_g1  463.00
    ECEF_Abd_TRINITY_DN10037_c0_g1  425.00
    ECEF_Abd_TRINITY_DN10037_c0_g2  554.00
    ECEF_Abd_TRINITY_DN10038_c0_g1  1041.00
    ECEF_Abd_TRINITY_DN10044_c0_g1  1605.00
    ECEF_Abd_TRINITY_DN10046_c0_g1  1211.00
    ECEF_Abd_TRINITY_DN10056_c0_g1  1544.00

With that preliminary input information, we are ready to examine GO term
enrichment with goseq\!

``` r
options(stringsAsFactors = FALSE) # this just always seems to help me. 
library("goseq")
```

    ## Loading required package: BiasedUrn

    ## Loading required package: geneLenDataBase

    ## 

``` r
library("dplyr")
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

I’m going to read in my DE analysis
results…

``` r
resdata <- read.table("C:/Users/cruth/Google Drive/Treehoppers/ResearchFiles/RNASeq/GeneExpression_2018/EC_DESeq2_resdata.tab", header=TRUE, sep="\t")
```

Not included in this walkthrough are the creation of the DESeq object,
doing the DE analysis, the creation of the resdata dataframe, and
writing that dataframe to a file. We only needed to do that once, and
it’s time intensive to repeat. I’ll just read in the data we saved
earlier. Right now, I’m only using it as a way to get the length of a
vector and to pick some example gene
ids.

### Here are all of the steps necessary for getting GO enrichment from a list of DE features.

1)  Make a vector the same length as the number of ids. All values will
    be 0.

<!-- end list -->

``` r
  gene.data <- integer(length=length(resdata$EC_id))
```

2)  Name the items in the vector with gene ids.

<!-- end list -->

``` r
names(gene.data) <- 
  resdata$EC_id
```

3)  Get a list of gene ids from SOMEWHERE – here, the top 30 values for
    ECEF Abd in resdata

<!-- end list -->

``` r
 topAbd <- arrange(resdata, ECEF_Abd) %>% top_n(30, ECEF_Abd)
 topAbd_ids <- topAbd$EC_id

AbdUp <- topAbd_ids
```

4)  Use a for loop to find which element in the genes vector has the
    same name as one of the gene ids in our list, and sets its value to
    1.

<!-- end list -->

``` r
for (item in AbdUp) {
  i <- which(names(gene.data)==item)
  gene.data[i] <- 1
}
```

5)  Check to make sure that the number of IDs you selected equals the
    number that are now set to 1 in genes

<!-- end list -->

``` r
table(gene.data)
```

    ## gene.data
    ##     0     1 
    ## 18645    30

``` r
length(AbdUp) == table(gene.data)[2]
```

    ##    1 
    ## TRUE

6)  Then omit NAs (probably skippable because of how we created the
    vector, but could be a problem for other means.)

<!-- end list -->

``` r
genes.list <- na.omit(gene.data)
```

7)  Get the feature length data and the category mapping data. Filter
    with dplyr to just be for the genes in our genes
list.

<!-- end list -->

``` r
feature.lengths <- read.table("../GOAnalysis/ECEF_Genes_Lengths.txt", header=TRUE, sep="\t")
cat.map <- read.table("../GOAnalysis/ECEF_GoTerm_Extended.txt", header=FALSE, sep="\t")
GoCats <- cat.map[,c(1,3)]

genes.length.data <- filter(feature.lengths, gene_id %in% names(genes.list))
```

8)  Make the length data into a vector with names.

<!-- end list -->

``` r
genes.bias.data <- genes.length.data$length
names(genes.bias.data) <- genes.length.data$gene_id
```

9)  Fit the probability weighting function and then plot
    it.

<!-- end list -->

``` r
pwf <- nullp(genes.list, bias.data = genes.bias.data, plot.fit = FALSE)
```

    ## Warning in pcls(G): initial point very close to some inequality constraints

``` r
plotPWF(pwf = pwf, binsize = 100) # <-- you can change the binsize if you like, bins of 200 genes seems good to me. 
```

![](GoSeq_Walkthrough_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

10) Perform the Wallenius test (use the goseq function) to get the most
    enriched go
cats

<!-- end list -->

``` r
GO.wall <- goseq(pwf, gene2cat = GoCats, method = "Wallenius", use_genes_without_cat = FALSE)
```

    ## Using manually entered categories.

    ## For 9464 genes, we could not find any categories. These genes will be excluded.

    ## To force their use, please run with use_genes_without_cat=TRUE (see documentation).

    ## This was the default behavior for version 1.15.1 and earlier.

    ## Calculating the p-values...

    ## 'select()' returned 1:1 mapping between keys and columns

Done. Save the data to a tab file to use for making graphics outputs.

``` r
write.table(GO.wall, "AbdUp_GOterms.wall.tab", sep="\t")
getwd()
```

    ## [1] "C:/Users/cruth/Google Drive/Treehoppers/ResearchFiles/RNASeq/GeneExpression_2018"

### And let’s take a look at the top twenty, shall we?

``` r
library(dplyr)
sorted.GO <- arrange(GO.wall, GO.wall$over_represented_pvalue) 

sorted.GO
```

    ##       category over_represented_pvalue under_represented_pvalue numDEInCat
    ## 1   GO:0042302            4.189228e-06               0.99999993          4
    ## 2   GO:0044421            9.433435e-05               0.99999450          5
    ## 3   GO:0042221            2.952566e-03               0.99949435          7
    ## 4   GO:0009605            3.455863e-03               0.99945735          6
    ## 5   GO:0043234            4.533398e-03               0.99905088          9
    ## 6   GO:0060361            2.519673e-02               0.99971993          1
    ## 7   GO:0008307            4.553681e-02               0.99904728          1
    ## 8   GO:0044700            7.039349e-02               0.97762676          6
    ## 9   GO:0044767            7.351987e-02               0.97361227          9
    ## 10  GO:0016043            7.457992e-02               0.97361351          8
    ## 11  GO:0007154            8.117343e-02               0.97317955          6
    ## 12  GO:0043228            1.258734e-01               0.95615901          5
    ## 13  GO:0048856            1.537212e-01               0.93388265          9
    ## 14  GO:0006457            1.579737e-01               0.98767105          1
    ## 15  GO:0009056            2.573098e-01               0.90736700          3
    ## 16  GO:0044707            2.870981e-01               0.85210019          9
    ## 17  GO:0044763            3.410801e-01               0.81794016         11
    ## 18  GO:0065008            3.751192e-01               0.83506550          3
    ## 19  GO:0044464            4.322528e-01               0.75551850         12
    ## 20  GO:0051641            4.872426e-01               0.78707330          2
    ## 21  GO:0003008            5.075241e-01               0.77086364          2
    ## 22  GO:0051674            5.115800e-01               0.84621710          1
    ## 23  GO:0044422            5.572335e-01               0.65063172          5
    ## 24  GO:0016787            5.880484e-01               0.65870884          3
    ## 25  GO:1901363            5.982722e-01               0.68979973          2
    ## 26  GO:0043167            6.025367e-01               0.77466352          1
    ## 27  GO:0097159            6.058126e-01               0.68239047          2
    ## 28  GO:0044710            6.310898e-01               0.61488192          3
    ## 29  GO:0051234            6.464740e-01               0.59844174          3
    ## 30  GO:0022892            6.647674e-01               0.71400020          1
    ## 31  GO:0022857            6.795716e-01               0.69793895          1
    ## 32  GO:0005515            6.998108e-01               0.53800614          3
    ## 33  GO:0044708            7.131659e-01               0.65887987          1
    ## 34  GO:0009058            7.615011e-01               0.50125295          2
    ## 35  GO:0044085            7.666233e-01               0.58825677          1
    ## 36  GO:0003006            7.730677e-01               0.57894812          1
    ## 37  GO:0006807            7.823141e-01               0.43236511          3
    ## 38  GO:0019953            8.274201e-01               0.49231479          1
    ## 39  GO:0043233            8.474156e-01               0.45620702          1
    ## 40  GO:0032504            8.711551e-01               0.40971361          1
    ## 41  GO:0016740            8.807352e-01               0.38968086          1
    ## 42  GO:0022414            8.994008e-01               0.34821617          1
    ## 43  GO:0006950            9.229466e-01               0.29044447          1
    ## 44  GO:0044238            9.572590e-01               0.11803186          4
    ## 45  GO:0044237            9.657626e-01               0.09874558          4
    ## 46  GO:0071704            9.683765e-01               0.09260354          4
    ## 47  GO:0043227            9.875353e-01               0.04691985          3
    ## 48  GO:0050789            9.985874e-01               0.00884843          2
    ## 49  GO:0097367            1.000000e+00               0.64544178          0
    ## 50  GO:0000150            1.000000e+00               0.99459436          0
    ## 51  GO:0000313            1.000000e+00               0.85746174          0
    ## 52  GO:0000989            1.000000e+00               0.76133692          0
    ## 53  GO:0000990            1.000000e+00               0.99639323          0
    ## 54  GO:0001776            1.000000e+00               0.99819509          0
    ## 55  GO:0001871            1.000000e+00               0.98904919          0
    ## 56  GO:0001909            1.000000e+00               0.99263586          0
    ## 57  GO:0001913            1.000000e+00               0.99459440          0
    ## 58  GO:0002209            1.000000e+00               0.99279544          0
    ## 59  GO:0002252            1.000000e+00               0.83771696          0
    ## 60  GO:0002262            1.000000e+00               0.96696752          0
    ## 61  GO:0002440            1.000000e+00               0.99096818          0
    ## 62  GO:0003700            1.000000e+00               0.36215194          0
    ## 63  GO:0003823            1.000000e+00               0.99076302          0
    ## 64  GO:0004129            1.000000e+00               0.98722873          0
    ## 65  GO:0004133            1.000000e+00               0.99638805          0
    ## 66  GO:0004362            1.000000e+00               0.99819509          0
    ## 67  GO:0004601            1.000000e+00               0.93451550          0
    ## 68  GO:0004708            1.000000e+00               0.99100277          0
    ## 69  GO:0004784            1.000000e+00               0.99099529          0
    ## 70  GO:0004791            1.000000e+00               0.99459349          0
    ## 71  GO:0004803            1.000000e+00               0.99100321          0
    ## 72  GO:0004872            1.000000e+00               0.40719251          0
    ## 73  GO:0004879            1.000000e+00               0.94418099          0
    ## 74  GO:0005057            1.000000e+00               0.87571020          0
    ## 75  GO:0005085            1.000000e+00               0.79668563          0
    ## 76  GO:0005200            1.000000e+00               0.91731711          0
    ## 77  GO:0005201            1.000000e+00               0.96605807          0
    ## 78  GO:0005212            1.000000e+00               0.99636539          0
    ## 79  GO:0005326            1.000000e+00               0.97052855          0
    ## 80  GO:0005487            1.000000e+00               0.97360397          0
    ## 81  GO:0005549            1.000000e+00               0.98978789          0
    ## 82  GO:0005911            1.000000e+00               0.81337077          0
    ## 83  GO:0006914            1.000000e+00               0.90109265          0
    ## 84  GO:0006955            1.000000e+00               0.47498234          0
    ## 85  GO:0007155            1.000000e+00               0.64596836          0
    ## 86  GO:0007622            1.000000e+00               0.71458229          0
    ## 87  GO:0007623            1.000000e+00               0.69155668          0
    ## 88  GO:0007624            1.000000e+00               0.99637325          0
    ## 89  GO:0007626            1.000000e+00               0.59477839          0
    ## 90  GO:0007631            1.000000e+00               0.92609550          0
    ## 91  GO:0008086            1.000000e+00               0.99639304          0
    ## 92  GO:0008144            1.000000e+00               0.95106680          0
    ## 93  GO:0008265            1.000000e+00               0.99459376          0
    ## 94  GO:0008283            1.000000e+00               0.63566350          0
    ## 95  GO:0008289            1.000000e+00               0.72357134          0
    ## 96  GO:0008384            1.000000e+00               0.99639324          0
    ## 97  GO:0008430            1.000000e+00               0.99639324          0
    ## 98  GO:0009607            1.000000e+00               0.39750159          0
    ## 99  GO:0009628            1.000000e+00               0.39997020          0
    ## 100 GO:0009719            1.000000e+00               0.64036754          0
    ## 101 GO:0009975            1.000000e+00               0.95376342          0
    ## 102 GO:0010166            1.000000e+00               0.99083807          0
    ## 103 GO:0010496            1.000000e+00               0.99261489          0
    ## 104 GO:0014823            1.000000e+00               0.99193481          0
    ## 105 GO:0016049            1.000000e+00               0.80143770          0
    ## 106 GO:0016247            1.000000e+00               0.97573628          0
    ## 107 GO:0016490            1.000000e+00               0.98913081          0
    ## 108 GO:0016491            1.000000e+00               0.44060250          0
    ## 109 GO:0016531            1.000000e+00               0.99639272          0
    ## 110 GO:0016533            1.000000e+00               0.99633008          0
    ## 111 GO:0016829            1.000000e+00               0.80057690          0
    ## 112 GO:0016853            1.000000e+00               0.85044328          0
    ## 113 GO:0017053            1.000000e+00               0.94084120          0
    ## 114 GO:0017056            1.000000e+00               0.97896389          0
    ## 115 GO:0019098            1.000000e+00               0.70285652          0
    ## 116 GO:0019239            1.000000e+00               0.97049897          0
    ## 117 GO:0019534            1.000000e+00               0.99638344          0
    ## 118 GO:0019825            1.000000e+00               0.99615042          0
    ## 119 GO:0019867            1.000000e+00               0.87910087          0
    ## 120 GO:0019882            1.000000e+00               0.98209792          0
    ## 121 GO:0019954            1.000000e+00               0.98741278          0
    ## 122 GO:0022403            1.000000e+00               0.85550412          0
    ## 123 GO:0022608            1.000000e+00               0.99819509          0
    ## 124 GO:0022611            1.000000e+00               0.99058191          0
    ## 125 GO:0030055            1.000000e+00               0.87923179          0
    ## 126 GO:0030234            1.000000e+00               0.50276907          0
    ## 127 GO:0030246            1.000000e+00               0.93780581          0
    ## 128 GO:0030297            1.000000e+00               0.99816357          0
    ## 129 GO:0030371            1.000000e+00               0.96889378          0
    ## 130 GO:0030545            1.000000e+00               0.98729284          0
    ## 131 GO:0031012            1.000000e+00               0.85933303          0
    ## 132 GO:0031594            1.000000e+00               0.86643442          0
    ## 133 GO:0032259            1.000000e+00               0.79959653          0
    ## 134 GO:0032451            1.000000e+00               0.97852572          0
    ## 135 GO:0032505            1.000000e+00               0.99276968          0
    ## 136 GO:0032947            1.000000e+00               0.99021627          0
    ## 137 GO:0032992            1.000000e+00               0.97999368          0
    ## 138 GO:0032993            1.000000e+00               0.93060483          0
    ## 139 GO:0032994            1.000000e+00               0.99456987          0
    ## 140 GO:0033036            1.000000e+00               0.26762745          0
    ## 141 GO:0033218            1.000000e+00               0.90530049          0
    ## 142 GO:0034518            1.000000e+00               0.99256174          0
    ## 143 GO:0034986            1.000000e+00               0.99819509          0
    ## 144 GO:0035172            1.000000e+00               0.97354795          0
    ## 145 GO:0035821            1.000000e+00               0.99639295          0
    ## 146 GO:0036094            1.000000e+00               0.64065460          0
    ## 147 GO:0036125            1.000000e+00               0.99639324          0
    ## 148 GO:0038023            1.000000e+00               0.48743446          0
    ## 149 GO:0042165            1.000000e+00               0.99242501          0
    ## 150 GO:0042267            1.000000e+00               0.99263586          0
    ## 151 GO:0042386            1.000000e+00               0.94699797          0
    ## 152 GO:0042445            1.000000e+00               0.85944048          0
    ## 153 GO:0042562            1.000000e+00               0.95030344          0
    ## 154 GO:0042910            1.000000e+00               0.94461077          0
    ## 155 GO:0042927            1.000000e+00               0.99639240          0
    ## 156 GO:0042954            1.000000e+00               0.99819323          0
    ## 157 GO:0043083            1.000000e+00               0.99819379          0
    ## 158 GO:0043176            1.000000e+00               0.99639209          0
    ## 159 GO:0043230            1.000000e+00               0.99459442          0
    ## 160 GO:0043235            1.000000e+00               0.91209980          0
    ## 161 GO:0043335            1.000000e+00               0.98381496          0
    ## 162 GO:0043473            1.000000e+00               0.83125131          0
    ## 163 GO:0043500            1.000000e+00               0.99819505          0
    ## 164 GO:0043627            1.000000e+00               0.97999121          0
    ## 165 GO:0044419            1.000000e+00               0.80406291          0
    ## 166 GO:0044425            1.000000e+00               0.07512958          0
    ## 167 GO:0044456            1.000000e+00               0.70544409          0
    ## 168 GO:0044706            1.000000e+00               0.79192367          0
    ## 169 GO:0044764            1.000000e+00               0.96767253          0
    ## 170 GO:0045153            1.000000e+00               0.99819509          0
    ## 171 GO:0045155            1.000000e+00               0.99819259          0
    ## 172 GO:0045174            1.000000e+00               0.99459445          0
    ## 173 GO:0045321            1.000000e+00               0.96385109          0
    ## 174 GO:0046790            1.000000e+00               0.99819502          0
    ## 175 GO:0048037            1.000000e+00               0.93612421          0
    ## 176 GO:0048475            1.000000e+00               0.95367497          0
    ## 177 GO:0048589            1.000000e+00               0.64880510          0
    ## 178 GO:0050840            1.000000e+00               0.99639324          0
    ## 179 GO:0050997            1.000000e+00               0.99053993          0
    ## 180 GO:0051183            1.000000e+00               0.98204018          0
    ## 181 GO:0051184            1.000000e+00               0.98192535          0
    ## 182 GO:0051235            1.000000e+00               0.83849188          0
    ## 183 GO:0051540            1.000000e+00               0.99255382          0
    ## 184 GO:0051606            1.000000e+00               0.65793252          0
    ## 185 GO:0051703            1.000000e+00               0.99421354          0
    ## 186 GO:0051705            1.000000e+00               0.67952246          0
    ## 187 GO:0051716            1.000000e+00               0.03537207          0
    ## 188 GO:0051920            1.000000e+00               0.98742165          0
    ## 189 GO:0060077            1.000000e+00               0.99459187          0
    ## 190 GO:0060090            1.000000e+00               0.93985084          0
    ## 191 GO:0065009            1.000000e+00               0.31133379          0
    ## 192 GO:0070161            1.000000e+00               0.77529132          0
    ## 193 GO:0070283            1.000000e+00               0.99622087          0
    ## 194 GO:0071554            1.000000e+00               0.99279563          0
    ## 195 GO:0072341            1.000000e+00               0.98172559          0
    ## 196 GO:0090079            1.000000e+00               0.96917980          0
    ## 197 GO:0090484            1.000000e+00               0.93587258          0
    ## 198 GO:0097060            1.000000e+00               0.91679524          0
    ## 199 GO:1901476            1.000000e+00               0.94122718          0
    ## 200 GO:1901505            1.000000e+00               0.94670767          0
    ## 201 GO:1901681            1.000000e+00               0.93384038          0
    ## 202 GO:0016874            1.000000e+00               0.62944884          0
    ## 203 GO:0003735            1.000000e+00               0.72673738          0
    ##     numInCat
    ## 1         62
    ## 2        269
    ## 3       1152
    ## 4        870
    ## 5       1995
    ## 6         14
    ## 7         25
    ## 8       1665
    ## 9       3053
    ## 10      2580
    ## 11      1725
    ## 12      1487
    ## 13      3496
    ## 14        92
    ## 15       959
    ## 16      3987
    ## 17      5224
    ## 18      1183
    ## 19      6019
    ## 20       872
    ## 21       903
    ## 22       380
    ## 23      2621
    ## 24      1607
    ## 25      1058
    ## 26       487
    ## 27      1072
    ## 28      1699
    ## 29      1741
    ## 30       579
    ## 31       601
    ## 32      1867
    ## 33       654
    ## 34      1414
    ## 35       754
    ## 36       768
    ## 37      2107
    ## 38       902
    ## 39       963
    ## 40      1045
    ## 41      1079
    ## 42      1161
    ## 43      1286
    ## 44      3725
    ## 45      3841
    ## 46      3882
    ## 47      3677
    ## 48      3823
    ## 49       236
    ## 50         3
    ## 51        83
    ## 52       146
    ## 53         2
    ## 54         1
    ## 55         6
    ## 56         4
    ## 57         3
    ## 58         4
    ## 59        96
    ## 60        18
    ## 61         5
    ## 62       534
    ## 63         5
    ## 64         7
    ## 65         2
    ## 66         1
    ## 67        37
    ## 68         5
    ## 69         5
    ## 70         3
    ## 71         5
    ## 72       476
    ## 73        31
    ## 74        72
    ## 75       121
    ## 76        47
    ## 77        19
    ## 78         2
    ## 79        16
    ## 80        14
    ## 81         5
    ## 82       110
    ## 83        56
    ## 84       396
    ## 85       233
    ## 86       181
    ## 87       198
    ## 88         2
    ## 89       277
    ## 90        41
    ## 91         2
    ## 92        27
    ## 93         3
    ## 94       240
    ## 95       173
    ## 96         2
    ## 97         2
    ## 98       485
    ## 99       481
    ## 100      237
    ## 101       26
    ## 102        5
    ## 103        4
    ## 104        4
    ## 105      119
    ## 106       13
    ## 107        6
    ## 108      434
    ## 109        2
    ## 110        2
    ## 111      121
    ## 112       86
    ## 113       33
    ## 114       11
    ## 115      189
    ## 116       16
    ## 117        2
    ## 118        2
    ## 119       70
    ## 120        9
    ## 121        7
    ## 122       83
    ## 123        1
    ## 124        5
    ## 125       69
    ## 126      363
    ## 127       35
    ## 128        1
    ## 129       17
    ## 130        7
    ## 131       82
    ## 132       77
    ## 133      121
    ## 134       12
    ## 135        4
    ## 136        5
    ## 137       11
    ## 138       39
    ## 139        3
    ## 140      689
    ## 141       54
    ## 142        4
    ## 143        1
    ## 144       14
    ## 145        2
    ## 146      239
    ## 147        2
    ## 148      383
    ## 149        4
    ## 150        4
    ## 151       28
    ## 152       80
    ## 153       27
    ## 154       31
    ## 155        2
    ## 156        1
    ## 157        1
    ## 158        2
    ## 159        3
    ## 160       50
    ## 161        9
    ## 162      100
    ## 163        1
    ## 164       11
    ## 165      118
    ## 166     1299
    ## 167      187
    ## 168      126
    ## 169       18
    ## 170        1
    ## 171        1
    ## 172        3
    ## 173       20
    ## 174        1
    ## 175       36
    ## 176       26
    ## 177      232
    ## 178        2
    ## 179        5
    ## 180       10
    ## 181       10
    ## 182       96
    ## 183        4
    ## 184      225
    ## 185        3
    ## 186      207
    ## 187     1641
    ## 188        7
    ## 189        3
    ## 190       32
    ## 191      609
    ## 192      135
    ## 193        2
    ## 194        4
    ## 195       10
    ## 196       17
    ## 197       36
    ## 198       47
    ## 199       33
    ## 200       30
    ## 201       37
    ## 202      248
    ## 203      171
    ##                                                                                                                                 term
    ## 1                                                                                                  structural constituent of cuticle
    ## 2                                                                                                          extracellular region part
    ## 3                                                                                                               response to chemical
    ## 4                                                                                                      response to external stimulus
    ## 5                                                                                                                               <NA>
    ## 6                                                                                                                             flight
    ## 7                                                                                                   structural constituent of muscle
    ## 8                                                                                                                               <NA>
    ## 9                                                                                                                               <NA>
    ## 10                                                                                                   cellular component organization
    ## 11                                                                                                                cell communication
    ## 12                                                                                                    non-membrane-bounded organelle
    ## 13                                                                                                  anatomical structure development
    ## 14                                                                                                                   protein folding
    ## 15                                                                                                                 catabolic process
    ## 16                                                                                                                              <NA>
    ## 17                                                                                                                              <NA>
    ## 18                                                                                                  regulation of biological quality
    ## 19                                                                                                                         cell part
    ## 20                                                                                                             cellular localization
    ## 21                                                                                                                    system process
    ## 22                                                                                                              localization of cell
    ## 23                                                                                                                    organelle part
    ## 24                                                                                                                hydrolase activity
    ## 25                                                                                                     heterocyclic compound binding
    ## 26                                                                                                                       ion binding
    ## 27                                                                                                   organic cyclic compound binding
    ## 28                                                                                                                              <NA>
    ## 29                                                                                                     establishment of localization
    ## 30                                                                                                                              <NA>
    ## 31                                                                                                transmembrane transporter activity
    ## 32                                                                                                                   protein binding
    ## 33                                                                                                                              <NA>
    ## 34                                                                                                              biosynthetic process
    ## 35                                                                                                     cellular component biogenesis
    ## 36                                                                                    developmental process involved in reproduction
    ## 37                                                                                               nitrogen compound metabolic process
    ## 38                                                                                                               sexual reproduction
    ## 39                                                                                                                   organelle lumen
    ## 40                                                                                               multicellular organism reproduction
    ## 41                                                                                                              transferase activity
    ## 42                                                                                                              reproductive process
    ## 43                                                                                                                response to stress
    ## 44                                                                                                         primary metabolic process
    ## 45                                                                                                        cellular metabolic process
    ## 46                                                                                               organic substance metabolic process
    ## 47                                                                                                        membrane-bounded organelle
    ## 48                                                                                                  regulation of biological process
    ## 49                                                                                                   carbohydrate derivative binding
    ## 50                                                                                                              recombinase activity
    ## 51                                                                                                               organellar ribosome
    ## 52                                                                                                                              <NA>
    ## 53                                                                                                                              <NA>
    ## 54                                                                                                             leukocyte homeostasis
    ## 55                                                                                                                   pattern binding
    ## 56                                                                                                   leukocyte mediated cytotoxicity
    ## 57                                                                                                      T cell mediated cytotoxicity
    ## 58                                                                                                       behavioral defense response
    ## 59                                                                                                           immune effector process
    ## 60                                                                                                          myeloid cell homeostasis
    ## 61                                                                               production of molecular mediator of immune response
    ## 62                                                                                         DNA-binding transcription factor activity
    ## 63                                                                                                                   antigen binding
    ## 64                                                                                                     cytochrome-c oxidase activity
    ## 65                                                                                              glycogen debranching enzyme activity
    ## 66                                                                                          glutathione-disulfide reductase activity
    ## 67                                                                                                               peroxidase activity
    ## 68                                                                                                        MAP kinase kinase activity
    ## 69                                                                                                     superoxide dismutase activity
    ## 70                                                                                          thioredoxin-disulfide reductase activity
    ## 71                                                                                                              transposase activity
    ## 72                                                                                                                              <NA>
    ## 73                                                                                                         nuclear receptor activity
    ## 74                                                                                                                              <NA>
    ## 75                                                                                        guanyl-nucleotide exchange factor activity
    ## 76                                                                                            structural constituent of cytoskeleton
    ## 77                                                                                       extracellular matrix structural constituent
    ## 78                                                                                                structural constituent of eye lens
    ## 79                                                                                             neurotransmitter transporter activity
    ## 80                                                                                                                              <NA>
    ## 81                                                                                                                   odorant binding
    ## 82                                                                                                                cell-cell junction
    ## 83                                                                                                                         autophagy
    ## 84                                                                                                                   immune response
    ## 85                                                                                                                     cell adhesion
    ## 86                                                                                                                 rhythmic behavior
    ## 87                                                                                                                  circadian rhythm
    ## 88                                                                                                                  ultradian rhythm
    ## 89                                                                                                               locomotory behavior
    ## 90                                                                                                                  feeding behavior
    ## 91                                                                            light-activated voltage-gated calcium channel activity
    ## 92                                                                                                                      drug binding
    ## 93                                                                                      Mo-molybdopterin cofactor sulfurase activity
    ## 94                                                                                                                cell proliferation
    ## 95                                                                                                                     lipid binding
    ## 96                                                                                                           IkappaB kinase activity
    ## 97                                                                                                                  selenium binding
    ## 98                                                                                                       response to biotic stimulus
    ## 99                                                                                                      response to abiotic stimulus
    ## 100                                                                                                  response to endogenous stimulus
    ## 101                                                                                                                 cyclase activity
    ## 102                                                                                                            wax metabolic process
    ## 103                                                                                                          intercellular transport
    ## 104                                                                                                             response to activity
    ## 105                                                                                                                      cell growth
    ## 106                                                                                                       channel regulator activity
    ## 107                                                                                   structural constituent of peritrophic membrane
    ## 108                                                                                                          oxidoreductase activity
    ## 109                                                                                                        copper chaperone activity
    ## 110                                                                                                         protein kinase 5 complex
    ## 111                                                                                                                   lyase activity
    ## 112                                                                                                               isomerase activity
    ## 113                                                                                                transcriptional repressor complex
    ## 114                                                                                           structural constituent of nuclear pore
    ## 115                                                                                                            reproductive behavior
    ## 116                                                                                                               deaminase activity
    ## 117                                                                                         toxin transmembrane transporter activity
    ## 118                                                                                                                   oxygen binding
    ## 119                                                                                                                   outer membrane
    ## 120                                                                                              antigen processing and presentation
    ## 121                                                                                                             asexual reproduction
    ## 122                                                                                                                 cell cycle phase
    ## 123                                                                                                  multicellular organism adhesion
    ## 124                                                                                                                 dormancy process
    ## 125                                                                                                          cell-substrate junction
    ## 126                                                                                                        enzyme regulator activity
    ## 127                                                                                                             carbohydrate binding
    ## 128                                                                transmembrane receptor protein tyrosine kinase activator activity
    ## 129                                                                                                   translation repressor activity
    ## 130                                                                                                      receptor regulator activity
    ## 131                                                                                                             extracellular matrix
    ## 132                                                                                                           neuromuscular junction
    ## 133                                                                                                                      methylation
    ## 134                                                                                                             demethylase activity
    ## 135                                                                                         reproduction of a single-celled organism
    ## 136                                                                                     protein-containing complex scaffold activity
    ## 137                                                                                                     protein-carbohydrate complex
    ## 138                                                                                                              protein-DNA complex
    ## 139                                                                                                            protein-lipid complex
    ## 140                                                                                                       macromolecule localization
    ## 141                                                                                                                    amide binding
    ## 142                                                                                                          RNA cap binding complex
    ## 143                                                                                                          iron chaperone activity
    ## 144                                                                                                           hemocyte proliferation
    ## 145                                                                       modification of morphology or physiology of other organism
    ## 146                                                                                                           small molecule binding
    ## 147                                                                                    fatty acid beta-oxidation multienzyme complex
    ## 148                                                                                                      signaling receptor activity
    ## 149                                                                                                         neurotransmitter binding
    ## 150                                                                                        natural killer cell mediated cytotoxicity
    ## 151                                                                                                         hemocyte differentiation
    ## 152                                                                                                        hormone metabolic process
    ## 153                                                                                                                  hormone binding
    ## 154                                                                                    xenobiotic transmembrane transporter activity
    ## 155                                                                                                                             <NA>
    ## 156                                                                                                 lipoprotein transporter activity
    ## 157                                                                                                                   synaptic cleft
    ## 158                                                                                                                    amine binding
    ## 159                                                                                                          extracellular organelle
    ## 160                                                                                                                 receptor complex
    ## 161                                                                                                                protein unfolding
    ## 162                                                                                                                     pigmentation
    ## 163                                                                                                                muscle adaptation
    ## 164                                                                                                             response to estrogen
    ## 165                                                                                       interspecies interaction between organisms
    ## 166                                                                                                                    membrane part
    ## 167                                                                                                                     synapse part
    ## 168                                                                                             multi-multicellular organism process
    ## 169                                                                                                  multi-organism cellular process
    ## 170                                electron transporter, transferring electrons within CoQH2-cytochrome c reductase complex activity
    ## 171 electron transporter, transferring electrons from CoQH2-cytochrome c reductase complex and cytochrome c oxidase complex activity
    ## 172                                                                                   glutathione dehydrogenase (ascorbate) activity
    ## 173                                                                                                             leukocyte activation
    ## 174                                                                                                                   virion binding
    ## 175                                                                                                                 cofactor binding
    ## 176                                                                                                                  coated membrane
    ## 177                                                                                                             developmental growth
    ## 178                                                                                                     extracellular matrix binding
    ## 179                                                                                                quaternary ammonium group binding
    ## 180                                                                                                                             <NA>
    ## 181                                                                                      cofactor transmembrane transporter activity
    ## 182                                                                                                          maintenance of location
    ## 183                                                                                                            metal cluster binding
    ## 184                                                                                                            detection of stimulus
    ## 185                                                                                       intraspecies interaction between organisms
    ## 186                                                                                                          multi-organism behavior
    ## 187                                                                                                    cellular response to stimulus
    ## 188                                                                                                           peroxiredoxin activity
    ## 189                                                                                                               inhibitory synapse
    ## 190                                                                                                       molecular adaptor activity
    ## 191                                                                                                 regulation of molecular function
    ## 192                                                                                                               anchoring junction
    ## 193                                                                                                                             <NA>
    ## 194                                                                                             cell wall organization or biogenesis
    ## 195                                                                                                      modified amino acid binding
    ## 196                                                                             translation regulator activity, nucleic acid binding
    ## 197                                                                                                                             <NA>
    ## 198                                                                                                                synaptic membrane
    ## 199                                                                                                                             <NA>
    ## 200                                                                       carbohydrate derivative transmembrane transporter activity
    ## 201                                                                                                          sulfur compound binding
    ## 202                                                                                                                  ligase activity
    ## 203                                                                                               structural constituent of ribosome
    ##     ontology
    ## 1         MF
    ## 2         CC
    ## 3         BP
    ## 4         BP
    ## 5       <NA>
    ## 6         BP
    ## 7         MF
    ## 8       <NA>
    ## 9       <NA>
    ## 10        BP
    ## 11        BP
    ## 12        CC
    ## 13        BP
    ## 14        BP
    ## 15        BP
    ## 16      <NA>
    ## 17      <NA>
    ## 18        BP
    ## 19        CC
    ## 20        BP
    ## 21        BP
    ## 22        BP
    ## 23        CC
    ## 24        MF
    ## 25        MF
    ## 26        MF
    ## 27        MF
    ## 28      <NA>
    ## 29        BP
    ## 30      <NA>
    ## 31        MF
    ## 32        MF
    ## 33      <NA>
    ## 34        BP
    ## 35        BP
    ## 36        BP
    ## 37        BP
    ## 38        BP
    ## 39        CC
    ## 40        BP
    ## 41        MF
    ## 42        BP
    ## 43        BP
    ## 44        BP
    ## 45        BP
    ## 46        BP
    ## 47        CC
    ## 48        BP
    ## 49        MF
    ## 50        MF
    ## 51        CC
    ## 52      <NA>
    ## 53      <NA>
    ## 54        BP
    ## 55        MF
    ## 56        BP
    ## 57        BP
    ## 58        BP
    ## 59        BP
    ## 60        BP
    ## 61        BP
    ## 62        MF
    ## 63        MF
    ## 64        MF
    ## 65        MF
    ## 66        MF
    ## 67        MF
    ## 68        MF
    ## 69        MF
    ## 70        MF
    ## 71        MF
    ## 72      <NA>
    ## 73        MF
    ## 74      <NA>
    ## 75        MF
    ## 76        MF
    ## 77        MF
    ## 78        MF
    ## 79        MF
    ## 80      <NA>
    ## 81        MF
    ## 82        CC
    ## 83        BP
    ## 84        BP
    ## 85        BP
    ## 86        BP
    ## 87        BP
    ## 88        BP
    ## 89        BP
    ## 90        BP
    ## 91        MF
    ## 92        MF
    ## 93        MF
    ## 94        BP
    ## 95        MF
    ## 96        MF
    ## 97        MF
    ## 98        BP
    ## 99        BP
    ## 100       BP
    ## 101       MF
    ## 102       BP
    ## 103       BP
    ## 104       BP
    ## 105       BP
    ## 106       MF
    ## 107       MF
    ## 108       MF
    ## 109       MF
    ## 110       CC
    ## 111       MF
    ## 112       MF
    ## 113       CC
    ## 114       MF
    ## 115       BP
    ## 116       MF
    ## 117       MF
    ## 118       MF
    ## 119       CC
    ## 120       BP
    ## 121       BP
    ## 122       BP
    ## 123       BP
    ## 124       BP
    ## 125       CC
    ## 126       MF
    ## 127       MF
    ## 128       MF
    ## 129       MF
    ## 130       MF
    ## 131       CC
    ## 132       CC
    ## 133       BP
    ## 134       MF
    ## 135       BP
    ## 136       MF
    ## 137       CC
    ## 138       CC
    ## 139       CC
    ## 140       BP
    ## 141       MF
    ## 142       CC
    ## 143       MF
    ## 144       BP
    ## 145       BP
    ## 146       MF
    ## 147       CC
    ## 148       MF
    ## 149       MF
    ## 150       BP
    ## 151       BP
    ## 152       BP
    ## 153       MF
    ## 154       MF
    ## 155     <NA>
    ## 156       MF
    ## 157       CC
    ## 158       MF
    ## 159       CC
    ## 160       CC
    ## 161       BP
    ## 162       BP
    ## 163       BP
    ## 164       BP
    ## 165       BP
    ## 166       CC
    ## 167       CC
    ## 168       BP
    ## 169       BP
    ## 170       MF
    ## 171       MF
    ## 172       MF
    ## 173       BP
    ## 174       MF
    ## 175       MF
    ## 176       CC
    ## 177       BP
    ## 178       MF
    ## 179       MF
    ## 180     <NA>
    ## 181       MF
    ## 182       BP
    ## 183       MF
    ## 184       BP
    ## 185       BP
    ## 186       BP
    ## 187       BP
    ## 188       MF
    ## 189       CC
    ## 190       MF
    ## 191       BP
    ## 192       CC
    ## 193     <NA>
    ## 194       BP
    ## 195       MF
    ## 196       MF
    ## 197     <NA>
    ## 198       CC
    ## 199     <NA>
    ## 200       MF
    ## 201       MF
    ## 202       MF
    ## 203       MF

``` r
#What are the top 20 over represented GOTerms?
sorted.GO[,c(1,6)]
```

    ##       category
    ## 1   GO:0042302
    ## 2   GO:0044421
    ## 3   GO:0042221
    ## 4   GO:0009605
    ## 5   GO:0043234
    ## 6   GO:0060361
    ## 7   GO:0008307
    ## 8   GO:0044700
    ## 9   GO:0044767
    ## 10  GO:0016043
    ## 11  GO:0007154
    ## 12  GO:0043228
    ## 13  GO:0048856
    ## 14  GO:0006457
    ## 15  GO:0009056
    ## 16  GO:0044707
    ## 17  GO:0044763
    ## 18  GO:0065008
    ## 19  GO:0044464
    ## 20  GO:0051641
    ## 21  GO:0003008
    ## 22  GO:0051674
    ## 23  GO:0044422
    ## 24  GO:0016787
    ## 25  GO:1901363
    ## 26  GO:0043167
    ## 27  GO:0097159
    ## 28  GO:0044710
    ## 29  GO:0051234
    ## 30  GO:0022892
    ## 31  GO:0022857
    ## 32  GO:0005515
    ## 33  GO:0044708
    ## 34  GO:0009058
    ## 35  GO:0044085
    ## 36  GO:0003006
    ## 37  GO:0006807
    ## 38  GO:0019953
    ## 39  GO:0043233
    ## 40  GO:0032504
    ## 41  GO:0016740
    ## 42  GO:0022414
    ## 43  GO:0006950
    ## 44  GO:0044238
    ## 45  GO:0044237
    ## 46  GO:0071704
    ## 47  GO:0043227
    ## 48  GO:0050789
    ## 49  GO:0097367
    ## 50  GO:0000150
    ## 51  GO:0000313
    ## 52  GO:0000989
    ## 53  GO:0000990
    ## 54  GO:0001776
    ## 55  GO:0001871
    ## 56  GO:0001909
    ## 57  GO:0001913
    ## 58  GO:0002209
    ## 59  GO:0002252
    ## 60  GO:0002262
    ## 61  GO:0002440
    ## 62  GO:0003700
    ## 63  GO:0003823
    ## 64  GO:0004129
    ## 65  GO:0004133
    ## 66  GO:0004362
    ## 67  GO:0004601
    ## 68  GO:0004708
    ## 69  GO:0004784
    ## 70  GO:0004791
    ## 71  GO:0004803
    ## 72  GO:0004872
    ## 73  GO:0004879
    ## 74  GO:0005057
    ## 75  GO:0005085
    ## 76  GO:0005200
    ## 77  GO:0005201
    ## 78  GO:0005212
    ## 79  GO:0005326
    ## 80  GO:0005487
    ## 81  GO:0005549
    ## 82  GO:0005911
    ## 83  GO:0006914
    ## 84  GO:0006955
    ## 85  GO:0007155
    ## 86  GO:0007622
    ## 87  GO:0007623
    ## 88  GO:0007624
    ## 89  GO:0007626
    ## 90  GO:0007631
    ## 91  GO:0008086
    ## 92  GO:0008144
    ## 93  GO:0008265
    ## 94  GO:0008283
    ## 95  GO:0008289
    ## 96  GO:0008384
    ## 97  GO:0008430
    ## 98  GO:0009607
    ## 99  GO:0009628
    ## 100 GO:0009719
    ## 101 GO:0009975
    ## 102 GO:0010166
    ## 103 GO:0010496
    ## 104 GO:0014823
    ## 105 GO:0016049
    ## 106 GO:0016247
    ## 107 GO:0016490
    ## 108 GO:0016491
    ## 109 GO:0016531
    ## 110 GO:0016533
    ## 111 GO:0016829
    ## 112 GO:0016853
    ## 113 GO:0017053
    ## 114 GO:0017056
    ## 115 GO:0019098
    ## 116 GO:0019239
    ## 117 GO:0019534
    ## 118 GO:0019825
    ## 119 GO:0019867
    ## 120 GO:0019882
    ## 121 GO:0019954
    ## 122 GO:0022403
    ## 123 GO:0022608
    ## 124 GO:0022611
    ## 125 GO:0030055
    ## 126 GO:0030234
    ## 127 GO:0030246
    ## 128 GO:0030297
    ## 129 GO:0030371
    ## 130 GO:0030545
    ## 131 GO:0031012
    ## 132 GO:0031594
    ## 133 GO:0032259
    ## 134 GO:0032451
    ## 135 GO:0032505
    ## 136 GO:0032947
    ## 137 GO:0032992
    ## 138 GO:0032993
    ## 139 GO:0032994
    ## 140 GO:0033036
    ## 141 GO:0033218
    ## 142 GO:0034518
    ## 143 GO:0034986
    ## 144 GO:0035172
    ## 145 GO:0035821
    ## 146 GO:0036094
    ## 147 GO:0036125
    ## 148 GO:0038023
    ## 149 GO:0042165
    ## 150 GO:0042267
    ## 151 GO:0042386
    ## 152 GO:0042445
    ## 153 GO:0042562
    ## 154 GO:0042910
    ## 155 GO:0042927
    ## 156 GO:0042954
    ## 157 GO:0043083
    ## 158 GO:0043176
    ## 159 GO:0043230
    ## 160 GO:0043235
    ## 161 GO:0043335
    ## 162 GO:0043473
    ## 163 GO:0043500
    ## 164 GO:0043627
    ## 165 GO:0044419
    ## 166 GO:0044425
    ## 167 GO:0044456
    ## 168 GO:0044706
    ## 169 GO:0044764
    ## 170 GO:0045153
    ## 171 GO:0045155
    ## 172 GO:0045174
    ## 173 GO:0045321
    ## 174 GO:0046790
    ## 175 GO:0048037
    ## 176 GO:0048475
    ## 177 GO:0048589
    ## 178 GO:0050840
    ## 179 GO:0050997
    ## 180 GO:0051183
    ## 181 GO:0051184
    ## 182 GO:0051235
    ## 183 GO:0051540
    ## 184 GO:0051606
    ## 185 GO:0051703
    ## 186 GO:0051705
    ## 187 GO:0051716
    ## 188 GO:0051920
    ## 189 GO:0060077
    ## 190 GO:0060090
    ## 191 GO:0065009
    ## 192 GO:0070161
    ## 193 GO:0070283
    ## 194 GO:0071554
    ## 195 GO:0072341
    ## 196 GO:0090079
    ## 197 GO:0090484
    ## 198 GO:0097060
    ## 199 GO:1901476
    ## 200 GO:1901505
    ## 201 GO:1901681
    ## 202 GO:0016874
    ## 203 GO:0003735
    ##                                                                                                                                 term
    ## 1                                                                                                  structural constituent of cuticle
    ## 2                                                                                                          extracellular region part
    ## 3                                                                                                               response to chemical
    ## 4                                                                                                      response to external stimulus
    ## 5                                                                                                                               <NA>
    ## 6                                                                                                                             flight
    ## 7                                                                                                   structural constituent of muscle
    ## 8                                                                                                                               <NA>
    ## 9                                                                                                                               <NA>
    ## 10                                                                                                   cellular component organization
    ## 11                                                                                                                cell communication
    ## 12                                                                                                    non-membrane-bounded organelle
    ## 13                                                                                                  anatomical structure development
    ## 14                                                                                                                   protein folding
    ## 15                                                                                                                 catabolic process
    ## 16                                                                                                                              <NA>
    ## 17                                                                                                                              <NA>
    ## 18                                                                                                  regulation of biological quality
    ## 19                                                                                                                         cell part
    ## 20                                                                                                             cellular localization
    ## 21                                                                                                                    system process
    ## 22                                                                                                              localization of cell
    ## 23                                                                                                                    organelle part
    ## 24                                                                                                                hydrolase activity
    ## 25                                                                                                     heterocyclic compound binding
    ## 26                                                                                                                       ion binding
    ## 27                                                                                                   organic cyclic compound binding
    ## 28                                                                                                                              <NA>
    ## 29                                                                                                     establishment of localization
    ## 30                                                                                                                              <NA>
    ## 31                                                                                                transmembrane transporter activity
    ## 32                                                                                                                   protein binding
    ## 33                                                                                                                              <NA>
    ## 34                                                                                                              biosynthetic process
    ## 35                                                                                                     cellular component biogenesis
    ## 36                                                                                    developmental process involved in reproduction
    ## 37                                                                                               nitrogen compound metabolic process
    ## 38                                                                                                               sexual reproduction
    ## 39                                                                                                                   organelle lumen
    ## 40                                                                                               multicellular organism reproduction
    ## 41                                                                                                              transferase activity
    ## 42                                                                                                              reproductive process
    ## 43                                                                                                                response to stress
    ## 44                                                                                                         primary metabolic process
    ## 45                                                                                                        cellular metabolic process
    ## 46                                                                                               organic substance metabolic process
    ## 47                                                                                                        membrane-bounded organelle
    ## 48                                                                                                  regulation of biological process
    ## 49                                                                                                   carbohydrate derivative binding
    ## 50                                                                                                              recombinase activity
    ## 51                                                                                                               organellar ribosome
    ## 52                                                                                                                              <NA>
    ## 53                                                                                                                              <NA>
    ## 54                                                                                                             leukocyte homeostasis
    ## 55                                                                                                                   pattern binding
    ## 56                                                                                                   leukocyte mediated cytotoxicity
    ## 57                                                                                                      T cell mediated cytotoxicity
    ## 58                                                                                                       behavioral defense response
    ## 59                                                                                                           immune effector process
    ## 60                                                                                                          myeloid cell homeostasis
    ## 61                                                                               production of molecular mediator of immune response
    ## 62                                                                                         DNA-binding transcription factor activity
    ## 63                                                                                                                   antigen binding
    ## 64                                                                                                     cytochrome-c oxidase activity
    ## 65                                                                                              glycogen debranching enzyme activity
    ## 66                                                                                          glutathione-disulfide reductase activity
    ## 67                                                                                                               peroxidase activity
    ## 68                                                                                                        MAP kinase kinase activity
    ## 69                                                                                                     superoxide dismutase activity
    ## 70                                                                                          thioredoxin-disulfide reductase activity
    ## 71                                                                                                              transposase activity
    ## 72                                                                                                                              <NA>
    ## 73                                                                                                         nuclear receptor activity
    ## 74                                                                                                                              <NA>
    ## 75                                                                                        guanyl-nucleotide exchange factor activity
    ## 76                                                                                            structural constituent of cytoskeleton
    ## 77                                                                                       extracellular matrix structural constituent
    ## 78                                                                                                structural constituent of eye lens
    ## 79                                                                                             neurotransmitter transporter activity
    ## 80                                                                                                                              <NA>
    ## 81                                                                                                                   odorant binding
    ## 82                                                                                                                cell-cell junction
    ## 83                                                                                                                         autophagy
    ## 84                                                                                                                   immune response
    ## 85                                                                                                                     cell adhesion
    ## 86                                                                                                                 rhythmic behavior
    ## 87                                                                                                                  circadian rhythm
    ## 88                                                                                                                  ultradian rhythm
    ## 89                                                                                                               locomotory behavior
    ## 90                                                                                                                  feeding behavior
    ## 91                                                                            light-activated voltage-gated calcium channel activity
    ## 92                                                                                                                      drug binding
    ## 93                                                                                      Mo-molybdopterin cofactor sulfurase activity
    ## 94                                                                                                                cell proliferation
    ## 95                                                                                                                     lipid binding
    ## 96                                                                                                           IkappaB kinase activity
    ## 97                                                                                                                  selenium binding
    ## 98                                                                                                       response to biotic stimulus
    ## 99                                                                                                      response to abiotic stimulus
    ## 100                                                                                                  response to endogenous stimulus
    ## 101                                                                                                                 cyclase activity
    ## 102                                                                                                            wax metabolic process
    ## 103                                                                                                          intercellular transport
    ## 104                                                                                                             response to activity
    ## 105                                                                                                                      cell growth
    ## 106                                                                                                       channel regulator activity
    ## 107                                                                                   structural constituent of peritrophic membrane
    ## 108                                                                                                          oxidoreductase activity
    ## 109                                                                                                        copper chaperone activity
    ## 110                                                                                                         protein kinase 5 complex
    ## 111                                                                                                                   lyase activity
    ## 112                                                                                                               isomerase activity
    ## 113                                                                                                transcriptional repressor complex
    ## 114                                                                                           structural constituent of nuclear pore
    ## 115                                                                                                            reproductive behavior
    ## 116                                                                                                               deaminase activity
    ## 117                                                                                         toxin transmembrane transporter activity
    ## 118                                                                                                                   oxygen binding
    ## 119                                                                                                                   outer membrane
    ## 120                                                                                              antigen processing and presentation
    ## 121                                                                                                             asexual reproduction
    ## 122                                                                                                                 cell cycle phase
    ## 123                                                                                                  multicellular organism adhesion
    ## 124                                                                                                                 dormancy process
    ## 125                                                                                                          cell-substrate junction
    ## 126                                                                                                        enzyme regulator activity
    ## 127                                                                                                             carbohydrate binding
    ## 128                                                                transmembrane receptor protein tyrosine kinase activator activity
    ## 129                                                                                                   translation repressor activity
    ## 130                                                                                                      receptor regulator activity
    ## 131                                                                                                             extracellular matrix
    ## 132                                                                                                           neuromuscular junction
    ## 133                                                                                                                      methylation
    ## 134                                                                                                             demethylase activity
    ## 135                                                                                         reproduction of a single-celled organism
    ## 136                                                                                     protein-containing complex scaffold activity
    ## 137                                                                                                     protein-carbohydrate complex
    ## 138                                                                                                              protein-DNA complex
    ## 139                                                                                                            protein-lipid complex
    ## 140                                                                                                       macromolecule localization
    ## 141                                                                                                                    amide binding
    ## 142                                                                                                          RNA cap binding complex
    ## 143                                                                                                          iron chaperone activity
    ## 144                                                                                                           hemocyte proliferation
    ## 145                                                                       modification of morphology or physiology of other organism
    ## 146                                                                                                           small molecule binding
    ## 147                                                                                    fatty acid beta-oxidation multienzyme complex
    ## 148                                                                                                      signaling receptor activity
    ## 149                                                                                                         neurotransmitter binding
    ## 150                                                                                        natural killer cell mediated cytotoxicity
    ## 151                                                                                                         hemocyte differentiation
    ## 152                                                                                                        hormone metabolic process
    ## 153                                                                                                                  hormone binding
    ## 154                                                                                    xenobiotic transmembrane transporter activity
    ## 155                                                                                                                             <NA>
    ## 156                                                                                                 lipoprotein transporter activity
    ## 157                                                                                                                   synaptic cleft
    ## 158                                                                                                                    amine binding
    ## 159                                                                                                          extracellular organelle
    ## 160                                                                                                                 receptor complex
    ## 161                                                                                                                protein unfolding
    ## 162                                                                                                                     pigmentation
    ## 163                                                                                                                muscle adaptation
    ## 164                                                                                                             response to estrogen
    ## 165                                                                                       interspecies interaction between organisms
    ## 166                                                                                                                    membrane part
    ## 167                                                                                                                     synapse part
    ## 168                                                                                             multi-multicellular organism process
    ## 169                                                                                                  multi-organism cellular process
    ## 170                                electron transporter, transferring electrons within CoQH2-cytochrome c reductase complex activity
    ## 171 electron transporter, transferring electrons from CoQH2-cytochrome c reductase complex and cytochrome c oxidase complex activity
    ## 172                                                                                   glutathione dehydrogenase (ascorbate) activity
    ## 173                                                                                                             leukocyte activation
    ## 174                                                                                                                   virion binding
    ## 175                                                                                                                 cofactor binding
    ## 176                                                                                                                  coated membrane
    ## 177                                                                                                             developmental growth
    ## 178                                                                                                     extracellular matrix binding
    ## 179                                                                                                quaternary ammonium group binding
    ## 180                                                                                                                             <NA>
    ## 181                                                                                      cofactor transmembrane transporter activity
    ## 182                                                                                                          maintenance of location
    ## 183                                                                                                            metal cluster binding
    ## 184                                                                                                            detection of stimulus
    ## 185                                                                                       intraspecies interaction between organisms
    ## 186                                                                                                          multi-organism behavior
    ## 187                                                                                                    cellular response to stimulus
    ## 188                                                                                                           peroxiredoxin activity
    ## 189                                                                                                               inhibitory synapse
    ## 190                                                                                                       molecular adaptor activity
    ## 191                                                                                                 regulation of molecular function
    ## 192                                                                                                               anchoring junction
    ## 193                                                                                                                             <NA>
    ## 194                                                                                             cell wall organization or biogenesis
    ## 195                                                                                                      modified amino acid binding
    ## 196                                                                             translation regulator activity, nucleic acid binding
    ## 197                                                                                                                             <NA>
    ## 198                                                                                                                synaptic membrane
    ## 199                                                                                                                             <NA>
    ## 200                                                                       carbohydrate derivative transmembrane transporter activity
    ## 201                                                                                                          sulfur compound binding
    ## 202                                                                                                                  ligase activity
    ## 203                                                                                               structural constituent of ribosome

Note – for some reason, many of the GOterms that EnTAP assigned are
deprecated, or maybe just not meant to be used the way we’re using them,
and so they show up as NAs. Those GOTerms can be searched on the Amigo
browser and you can find a redirect to the current ontology version. For
example, <GO:0043234> was replaced by <GO:0032991>, protein-containing
complex, a cellular component term.

That’s all\! The hard part is deciding *which* set of genes to examine.
