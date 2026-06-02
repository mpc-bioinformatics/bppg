# test aggregateReplicates

    Code
      D1
    Output
      class: SummarizedExperiment 
      dim: 87 9 
      metadata(0):
      assays(1): intensities
      rownames(87): AAADALSDLEIK AAADALSDLEIKDSK ... NMVSDIQEATK VLGIDAGEEK
      rowData names(1): Sequence
      colnames(9): 1 2 ... 8 9
      colData names(1): group

---

    Code
      SummarizedExperiment::assays(D1)$intensities
    Output
                                               1         2         3         4
      AAADALSDLEIK                            NA        NA        NA        NA
      AAADALSDLEIKDSK                    6983400  13878000   8056933   7299467
      AEWALR                             4981633   6038733   5616700   5051567
      DEGLHTDFACLLFAHLK                  2926297   4196733   5006267   4604733
      DIHDWNNR                                NA        NA   1805123        NA
      ELETLREENR                              NA   5833967        NA        NA
      ESEFLFNAIHTIPEIGEK                34531333  36904333  34725000  32236000
      GMMPGLTFSNELICR                    8221367   7740633   8593500   9115600
      IVTEAVEIEQR                       19285000  19356333  23024333  18153667
      LLVAFGNK                          12208000   7570000  13673333  10312200
      LLVAFGNKK                               NA   1930200        NA        NA
      NKPDPAIVEK                        16994333  19315333  18701000  15761667
      TNFFEK                             6472733   5585067   7092900   6198100
      TVLFPIK                           13587333  11388667  14516667  11367667
      VENPFDFMENISLAGK                  18441667  18036333  17716000  16306667
      WIQDADALFGER                      18621333  17592667  19626667  20071333
      YFLDALPVALLGMNADLMNQYVEFVADR      28552333  24190667  27130000  24436333
      AANLGGVAVSGLEMAQNSQK               7682633  10221667  10057200   9901233
      DAVWFGPPK                               NA        NA        NA        NA
      EIGYLFGAYR                              NA  25966000        NA        NA
      FHPSVNLSILK                             NA   5193033   7085767   6011367
      FLGFEQIFK                         40834000  49968333  39339333  40588000
      GANIASFVMVADAMLDQGDVF             56648667  85393000  50223333  54703333
      GCIISETGITSEQIHDIASAK              8447933   8636433   8286833   7580700
      GGLCVDLK                                NA        NA  12132300  10199767
      ICYAFMR                            9630700  12955667   8691367   9041467
      NSWEGVLTGK                        13002667        NA  15192000  12682667
      SLEEIVDEYSTFSESK                   9217967   9902000   9259467   9221533
      VLPIVSVPER                        38884333  36686667  47458667  41311000
      VTISGSGNVAQYAALK                   3394137        NA   2970800        NA
      VTWENDNGEQEVAQGYR                       NA        NA        NA   3155933
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK        NA        NA        NA        NA
      AANLGGVAVSGLEMAQNSQR              29756667  27488333  32588000  33017667
      ALVAQGVK                          12393667  14875333  17112000        NA
      FIAEGSNMGSTPEAIAVFETAR            36095000  44362333  42719333  36686000
      GANIASFIK                         17540333  19988667        NA  18060000
      GCIISETGITSEQVADISSAK              3639367   3074667   3396200   3512567
      HIGQDTDVPAGDIGVGGR                31470333  36311333  35960000  31697000
      IMINCFNECIDYAK                    13701000  13295000  12435333  13489667
      ITWTSER                           18395667  12407563  21792667        NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     63609333  76521333  89660000  89529333
      SLEQIVNEYSTFSENK                  60231667  50030333  58250667  54767667
      STATGPSEAVWYGPPK                  43254667  44978667  48126667  47938667
      VDIALPCATQNEVSGEEAK               48229333  52357667  55121667  54851000
      VIELGGTVVSLSDSK                   15023333  13853667  16089000  13941000
      VQYIAGARPWTHVQK                         NA        NA   4688400        NA
      VTWENDKGEQEVAQGYR                 35064000  45248000  34321333  33438333
      AAGLTAAYAR                        59266000  63375000  69987667  55438667
      APEAEQVLSAAATFPIAQPATDVEAR        29290000  22954667  33818333  27877000
      AVQDNGESAFR                       12264000  14980667  14452333  13343667
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     28282000  43794000  24112667  24053667
      GFTLAEVK                          29186333  34311667  30924667  30637000
      IAPRPLDLLRPVVR                    17974000  27753000  18007333  12850033
      IIVFPR                            49650333  31246033  58704333  55813333
      NQEIFDANVQR                       86158000  83463667  93019333  86492000
      TIGIAVDHR                         25202000  19952000  28531000  22916000
      VHFDQAGK                          10391167        NA  12552333        NA
      VHFDQAGKK                               NA        NA        NA        NA
      ANELLINVK                        172740000 174836667 194686667 151003333
      ANGTTVLVGMPAGAK                         NA        NA        NA        NA
      ATDGGAHGVINVSVSEAAIEASTR         113900000  78513000 101964000  84256667
      CCSDVFNQVVK                       49396000  59550333  52326000  52362333
      DIVGAVLK                                NA        NA  56520000        NA
      EALDFFAR                         188476667 206550000 208380000 195696667
      EKDIVGAVLK                        39416000  35036000  36228000  32234667
      GVIFYESHGK                        95035333 118343333 101534667  87530667
      IGDYAGIK                         136606667 140413333 174476667 135713333
      LPLVGGHEGAGVVVGMGENVK            128880000 119036000 119073333 123006667
      SANLMAGHWVAISGAAGGLGSLAVQYAK      27385000  29500333  38130000  29294000
      SIGGEVFIDFTK                      65106667  75924667  69992667  65335667
      SIPETQK                           27845000  42305000  21272333  28597333
      SISIVGSYVGNR                     266830000 275783333 297603333 259220000
      VLGIDGGEGK                              NA  11413667  28787333  20040000
      VLGIDGGEGKEELFR                  303350000 383026667 300430000 303976667
      VVGLSTLPEIYEK                    519833333 550590000 552153333 513773333
      YSGVCHTDLHAWHGDWPLPVK             98361000 107558667 109023333 146366667
      ANGTVVLVGLPAGAK                    4686733   5487200   5100067   4686033
      CSSDVFNHVVK                             NA   4663567        NA        NA
      DIPVPKPKPNELLINVK                       NA        NA        NA        NA
      VVGLSSLPEIYEK                     37969333  37831333  43313667  36479667
      DIPVPEPKPNEILINVK                 19119667  21145000  24534000  22568333
      EALDFFSR                           6554133   7199000   6608867   7492567
      GVIFYENK                           5322033        NA   5378500   5282367
      IQQGTDLAEVAPILCAGVTVYK             7426767   8251100   7680167   7360667
      IVGLSELPK                         17489667  20732000  19205667  19480667
      NMVSDIQEATK                        6085200        NA   6700833   6403567
      VLGIDAGEEK                              NA        NA        NA        NA
                                               5         6         7         8
      AAADALSDLEIK                            NA        NA        NA        NA
      AAADALSDLEIKDSK                   12349667   6444600        NA   6948833
      AEWALR                             5463800   5515467   5298200   5565500
      DEGLHTDFACLLFAHLK                       NA        NA        NA        NA
      DIHDWNNR                                NA        NA        NA        NA
      ELETLREENR                         4083867   5886800        NA   4740933
      ESEFLFNAIHTIPEIGEK                29790333  38798000  32100000  34559000
      GMMPGLTFSNELICR                         NA   8444400  10313933   9266667
      IVTEAVEIEQR                       18239667  26923333  17141667  18727000
      LLVAFGNK                           6906700  16992333        NA   8820667
      LLVAFGNKK                               NA        NA        NA        NA
      NKPDPAIVEK                        18437333  25011000  16730667  20630333
      TNFFEK                             5954433   6293500   6691533   6443500
      TVLFPIK                           10303733  16685333  11555667  11008000
      VENPFDFMENISLAGK                  17959667  14929000  18071000  17308333
      WIQDADALFGER                      15586667  16600667  20252667  16469667
      YFLDALPVALLGMNADLMNQYVEFVADR      23871333  27493667  23339667  26526333
      AANLGGVAVSGLEMAQNSQK               9044833  10327167   8229667  10649333
      DAVWFGPPK                               NA        NA        NA        NA
      EIGYLFGAYR                              NA        NA        NA        NA
      FHPSVNLSILK                        5620633   9841467   5110567   5119667
      FLGFEQIFK                         45476667  43832000  41277000  46683000
      GANIASFVMVADAMLDQGDVF             78859000  48051667  52853667  72607667
      GCIISETGITSEQIHDIASAK              8219633   9209100   7333567   8369633
      GGLCVDLK                          11244000  13672667        NA        NA
      ICYAFMR                           11167000  10049767   9564967  12442667
      NSWEGVLTGK                        13549000  16621000  11919667  14179667
      SLEEIVDEYSTFSESK                   8420900   6821100   8689267   9331000
      VLPIVSVPER                        35490333  46417667  44018333  43666333
      VTISGSGNVAQYAALK                        NA        NA        NA        NA
      VTWENDNGEQEVAQGYR                  2954267        NA   3551167   2916833
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK        NA        NA        NA        NA
      AANLGGVAVSGLEMAQNSQR              36991000  29288667  25432567  27723167
      ALVAQGVK                          13045000  16788667  12481000  14288667
      FIAEGSNMGSTPEAIAVFETAR            37530000  30513667  39268333  42755333
      GANIASFIK                         18407667  21705667  18207333  22172000
      GCIISETGITSEQVADISSAK              3784100   4597000   3756367   3777567
      HIGQDTDVPAGDIGVGGR                33837000  46534667  30756667  38132333
      IMINCFNECIDYAK                    13699667  14051667  13428333  14036000
      ITWTSER                           17077667        NA  19126000  19846000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR    100111333  36039333  84677667 100917333
      SLEQIVNEYSTFSENK                  55748667  58444333  58856333  54194667
      STATGPSEAVWYGPPK                  43205000        NA  47553667  47282333
      VDIALPCATQNEVSGEEAK               47385000  56077000  52440667  57681333
      VIELGGTVVSLSDSK                   13408000  15118333  13848667  15204000
      VQYIAGARPWTHVQK                    8435600        NA   9249400        NA
      VTWENDKGEQEVAQGYR                 39279000  42152000  35194000  41293000
      AAGLTAAYAR                        59486333  74013333  57707667  68245000
      APEAEQVLSAAATFPIAQPATDVEAR        17825667  32206667  27982333  17877000
      AVQDNGESAFR                        8764893  15664000  13687667  14547000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     27862333  15260667  29235333  27285000
      GFTLAEVK                          32141667  35082667  29205667  34466000
      IAPRPLDLLRPVVR                    20643333  16119667  15304533  28587333
      IIVFPR                            45084333  61537667  54826667  51289000
      NQEIFDANVQR                       83184333  96626333  93143000  93823667
      TIGIAVDHR                         21659667        NA  27198333  24122333
      VHFDQAGK                          10726667  24030333   9422000  11583333
      VHFDQAGKK                               NA        NA        NA        NA
      ANELLINVK                        163163333 242573333 164803333 168390000
      ANGTTVLVGMPAGAK                         NA  22699667  20349000        NA
      ATDGGAHGVINVSVSEAAIEASTR          73148667 141490000  90086667  94831333
      CCSDVFNQVVK                       54651667  57587667  49942000  53650000
      DIVGAVLK                                NA  62694000        NA  54120667
      EALDFFAR                         191590000 216690000 197646667 214413333
      EKDIVGAVLK                        35959667  46760667  33572333  37171333
      GVIFYESHGK                        96373333 138540000  89986000 100885000
      IGDYAGIK                                NA 169593333 138740000 131976667
      LPLVGGHEGAGVVVGMGENVK            108108667 151093333 128256667 141380000
      SANLMAGHWVAISGAAGGLGSLAVQYAK      28973667  27642667  21001333  36458000
      SIGGEVFIDFTK                      59080333  61164667  63580000  69802667
      SIPETQK                           25557333  12179000  28931667  22287333
      SISIVGSYVGNR                     250753333 335060000 273083333 285616667
      VLGIDGGEGK                        10551733  28280333  20313000   9879767
      VLGIDGGEGKEELFR                  327840000 362080000 302440000 357786667
      VVGLSTLPEIYEK                    513980000 554896667 540730000 570326667
      YSGVCHTDLHAWHGDWPLPVK            106637667 114550000 136323333 113767000
      ANGTVVLVGLPAGAK                    4931000   5919167   4421467   5358467
      CSSDVFNHVVK                             NA   4472800        NA        NA
      DIPVPKPKPNELLINVK                       NA        NA        NA        NA
      VVGLSSLPEIYEK                     34889000  40066667  39909000  35714000
      DIPVPEPKPNEILINVK                 14027917  22241333  20571667        NA
      EALDFFSR                                NA        NA   6708733   8074733
      GVIFYENK                                NA        NA        NA   5625633
      IQQGTDLAEVAPILCAGVTVYK             8412900   7482333   8204433   8673200
      IVGLSELPK                         18741000  22700333  19164000  21740667
      NMVSDIQEATK                        5964300   5878433   6147833   6105633
      VLGIDAGEEK                              NA        NA        NA        NA
                                               9
      AAADALSDLEIK                            NA
      AAADALSDLEIKDSK                         NA
      AEWALR                             5849200
      DEGLHTDFACLLFAHLK                  9592133
      DIHDWNNR                                NA
      ELETLREENR                         8287433
      ESEFLFNAIHTIPEIGEK                34887333
      GMMPGLTFSNELICR                         NA
      IVTEAVEIEQR                       26413667
      LLVAFGNK                           7729900
      LLVAFGNKK                          2905600
      NKPDPAIVEK                        21057000
      TNFFEK                             6567667
      TVLFPIK                           11048667
      VENPFDFMENISLAGK                  15313667
      WIQDADALFGER                      15994000
      YFLDALPVALLGMNADLMNQYVEFVADR      31549333
      AANLGGVAVSGLEMAQNSQK                    NA
      DAVWFGPPK                               NA
      EIGYLFGAYR                        28923000
      FHPSVNLSILK                       10438367
      FLGFEQIFK                         61947667
      GANIASFVMVADAMLDQGDVF            108181667
      GCIISETGITSEQIHDIASAK             10201567
      GGLCVDLK                          11996667
      ICYAFMR                           10883667
      NSWEGVLTGK                        14516333
      SLEEIVDEYSTFSESK                   8719967
      VLPIVSVPER                        35614667
      VTISGSGNVAQYAALK                        NA
      VTWENDNGEQEVAQGYR                  3404000
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK   6988267
      AANLGGVAVSGLEMAQNSQR              31534667
      ALVAQGVK                                NA
      FIAEGSNMGSTPEAIAVFETAR            43571333
      GANIASFIK                         18624667
      GCIISETGITSEQVADISSAK                   NA
      HIGQDTDVPAGDIGVGGR                41499333
      IMINCFNECIDYAK                    14511667
      ITWTSER                                 NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR           NA
      SLEQIVNEYSTFSENK                  37703667
      STATGPSEAVWYGPPK                  43159333
      VDIALPCATQNEVSGEEAK               35875467
      VIELGGTVVSLSDSK                   17355333
      VQYIAGARPWTHVQK                         NA
      VTWENDKGEQEVAQGYR                 57630333
      AAGLTAAYAR                        58873333
      APEAEQVLSAAATFPIAQPATDVEAR         3825933
      AVQDNGESAFR                       17361333
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     30293667
      GFTLAEVK                          36347333
      IAPRPLDLLRPVVR                    14247333
      IIVFPR                            16913867
      NQEIFDANVQR                       84064333
      TIGIAVDHR                         31816000
      VHFDQAGK                          15285000
      VHFDQAGKK                          2248033
      ANELLINVK                        182210000
      ANGTTVLVGMPAGAK                   21371667
      ATDGGAHGVINVSVSEAAIEASTR         105391667
      CCSDVFNQVVK                       69164000
      DIVGAVLK                          59953000
      EALDFFAR                         216606667
      EKDIVGAVLK                        43935667
      GVIFYESHGK                       209076667
      IGDYAGIK                                NA
      LPLVGGHEGAGVVVGMGENVK            124833333
      SANLMAGHWVAISGAAGGLGSLAVQYAK      50760333
      SIGGEVFIDFTK                      66369333
      SIPETQK                           41765333
      SISIVGSYVGNR                     278913333
      VLGIDGGEGK                        10213167
      VLGIDGGEGKEELFR                  467486667
      VVGLSTLPEIYEK                    549193333
      YSGVCHTDLHAWHGDWPLPVK            105074333
      ANGTVVLVGLPAGAK                    5926733
      CSSDVFNHVVK                        5493500
      DIPVPKPKPNELLINVK                  6552133
      VVGLSSLPEIYEK                     34009333
      DIPVPEPKPNEILINVK                 17775133
      EALDFFSR                           6977467
      GVIFYENK                           7224433
      IQQGTDLAEVAPILCAGVTVYK             8271300
      IVGLSELPK                         17076000
      NMVSDIQEATK                        7632267
      VLGIDAGEEK                              NA

---

    Code
      as.data.frame(tail(SummarizedExperiment::rowData(D1), n = 1000))
    Output
                                                               Sequence
      AAADALSDLEIK                                         AAADALSDLEIK
      AAADALSDLEIKDSK                                   AAADALSDLEIKDSK
      AEWALR                                                     AEWALR
      DEGLHTDFACLLFAHLK                               DEGLHTDFACLLFAHLK
      DIHDWNNR                                                 DIHDWNNR
      ELETLREENR                                             ELETLREENR
      ESEFLFNAIHTIPEIGEK                             ESEFLFNAIHTIPEIGEK
      GMMPGLTFSNELICR                                   GMMPGLTFSNELICR
      IVTEAVEIEQR                                           IVTEAVEIEQR
      LLVAFGNK                                                 LLVAFGNK
      LLVAFGNKK                                               LLVAFGNKK
      NKPDPAIVEK                                             NKPDPAIVEK
      TNFFEK                                                     TNFFEK
      TVLFPIK                                                   TVLFPIK
      VENPFDFMENISLAGK                                 VENPFDFMENISLAGK
      WIQDADALFGER                                         WIQDADALFGER
      YFLDALPVALLGMNADLMNQYVEFVADR         YFLDALPVALLGMNADLMNQYVEFVADR
      AANLGGVAVSGLEMAQNSQK                         AANLGGVAVSGLEMAQNSQK
      DAVWFGPPK                                               DAVWFGPPK
      EIGYLFGAYR                                             EIGYLFGAYR
      FHPSVNLSILK                                           FHPSVNLSILK
      FLGFEQIFK                                               FLGFEQIFK
      GANIASFVMVADAMLDQGDVF                       GANIASFVMVADAMLDQGDVF
      GCIISETGITSEQIHDIASAK                       GCIISETGITSEQIHDIASAK
      GGLCVDLK                                                 GGLCVDLK
      ICYAFMR                                                   ICYAFMR
      NSWEGVLTGK                                             NSWEGVLTGK
      SLEEIVDEYSTFSESK                                 SLEEIVDEYSTFSESK
      VLPIVSVPER                                             VLPIVSVPER
      VTISGSGNVAQYAALK                                 VTISGSGNVAQYAALK
      VTWENDNGEQEVAQGYR                               VTWENDNGEQEVAQGYR
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK YVAGARPWTHVSNVDIALPCATQNEVSGDEAK
      AANLGGVAVSGLEMAQNSQR                         AANLGGVAVSGLEMAQNSQR
      ALVAQGVK                                                 ALVAQGVK
      FIAEGSNMGSTPEAIAVFETAR                     FIAEGSNMGSTPEAIAVFETAR
      GANIASFIK                                               GANIASFIK
      GCIISETGITSEQVADISSAK                       GCIISETGITSEQVADISSAK
      HIGQDTDVPAGDIGVGGR                             HIGQDTDVPAGDIGVGGR
      IMINCFNECIDYAK                                     IMINCFNECIDYAK
      ITWTSER                                                   ITWTSER
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR       SEPEFQQAYEEVVSSLEDSTLFEQHPEYR
      SLEQIVNEYSTFSENK                                 SLEQIVNEYSTFSENK
      STATGPSEAVWYGPPK                                 STATGPSEAVWYGPPK
      VDIALPCATQNEVSGEEAK                           VDIALPCATQNEVSGEEAK
      VIELGGTVVSLSDSK                                   VIELGGTVVSLSDSK
      VQYIAGARPWTHVQK                                   VQYIAGARPWTHVQK
      VTWENDKGEQEVAQGYR                               VTWENDKGEQEVAQGYR
      AAGLTAAYAR                                             AAGLTAAYAR
      APEAEQVLSAAATFPIAQPATDVEAR             APEAEQVLSAAATFPIAQPATDVEAR
      AVQDNGESAFR                                           AVQDNGESAFR
      DGKAPEAEQVLSAAATFPIAQPATDVEAR       DGKAPEAEQVLSAAATFPIAQPATDVEAR
      GFTLAEVK                                                 GFTLAEVK
      IAPRPLDLLRPVVR                                     IAPRPLDLLRPVVR
      IIVFPR                                                     IIVFPR
      NQEIFDANVQR                                           NQEIFDANVQR
      TIGIAVDHR                                               TIGIAVDHR
      VHFDQAGK                                                 VHFDQAGK
      VHFDQAGKK                                               VHFDQAGKK
      ANELLINVK                                               ANELLINVK
      ANGTTVLVGMPAGAK                                   ANGTTVLVGMPAGAK
      ATDGGAHGVINVSVSEAAIEASTR                 ATDGGAHGVINVSVSEAAIEASTR
      CCSDVFNQVVK                                           CCSDVFNQVVK
      DIVGAVLK                                                 DIVGAVLK
      EALDFFAR                                                 EALDFFAR
      EKDIVGAVLK                                             EKDIVGAVLK
      GVIFYESHGK                                             GVIFYESHGK
      IGDYAGIK                                                 IGDYAGIK
      LPLVGGHEGAGVVVGMGENVK                       LPLVGGHEGAGVVVGMGENVK
      SANLMAGHWVAISGAAGGLGSLAVQYAK         SANLMAGHWVAISGAAGGLGSLAVQYAK
      SIGGEVFIDFTK                                         SIGGEVFIDFTK
      SIPETQK                                                   SIPETQK
      SISIVGSYVGNR                                         SISIVGSYVGNR
      VLGIDGGEGK                                             VLGIDGGEGK
      VLGIDGGEGKEELFR                                   VLGIDGGEGKEELFR
      VVGLSTLPEIYEK                                       VVGLSTLPEIYEK
      YSGVCHTDLHAWHGDWPLPVK                       YSGVCHTDLHAWHGDWPLPVK
      ANGTVVLVGLPAGAK                                   ANGTVVLVGLPAGAK
      CSSDVFNHVVK                                           CSSDVFNHVVK
      DIPVPKPKPNELLINVK                               DIPVPKPKPNELLINVK
      VVGLSSLPEIYEK                                       VVGLSSLPEIYEK
      DIPVPEPKPNEILINVK                               DIPVPEPKPNEILINVK
      EALDFFSR                                                 EALDFFSR
      GVIFYENK                                                 GVIFYENK
      IQQGTDLAEVAPILCAGVTVYK                     IQQGTDLAEVAPILCAGVTVYK
      IVGLSELPK                                               IVGLSELPK
      NMVSDIQEATK                                           NMVSDIQEATK
      VLGIDAGEEK                                             VLGIDAGEEK

---

    Code
      as.data.frame(tail(SummarizedExperiment::colData(D1), n = 1000))
    Output
        group
      1     1
      2     2
      3     3
      4     4
      5     5
      6     6
      7     7
      8     8
      9     9

---

    Code
      D2
    Output
      class: SummarizedExperiment 
      dim: 87 9 
      metadata(0):
      assays(1): intensities
      rownames(87): AAADALSDLEIK AAADALSDLEIKDSK ... NMVSDIQEATK VLGIDAGEEK
      rowData names(1): Sequence
      colnames(9): 1 2 ... 8 9
      colData names(1): group

---

    Code
      SummarizedExperiment::assays(D2)$intensities
    Output
                                               1         2         3         4
      AAADALSDLEIK                            NA        NA        NA        NA
      AAADALSDLEIKDSK                    7759800  13943000   7192900   7339300
      AEWALR                             4987700   6124700   6064200   4861600
      DEGLHTDFACLLFAHLK                  3640300   4267500   4885700   4938900
      DIHDWNNR                           1911100   2426550   2221900        NA
      ELETLREENR                              NA   5610400        NA        NA
      ESEFLFNAIHTIPEIGEK                34234000  36455000  33646000  35195000
      GMMPGLTFSNELICR                    8162900   7910400   8268200   9366600
      IVTEAVEIEQR                       19168000  19237000  24877000  17987000
      LLVAFGNK                          11947000   7773000  13451000  10335000
      LLVAFGNKK                               NA   1866500        NA        NA
      NKPDPAIVEK                        16866000  19391000  19010000  14737000
      TNFFEK                             7189600   5396000   7485700   6161400
      TVLFPIK                           13675000  10894000  15000000  11354000
      VENPFDFMENISLAGK                  18635000  18099000  17240000  16315000
      WIQDADALFGER                      18409000  17600000  19526000  19619000
      YFLDALPVALLGMNADLMNQYVEFVADR      28416000  23392000  26404000  24870000
      AANLGGVAVSGLEMAQNSQK               9531300  10051000   9545800  10264000
      DAVWFGPPK                               NA        NA        NA        NA
      EIGYLFGAYR                              NA  26668000        NA        NA
      FHPSVNLSILK                        4628050   5235400   7422800   5560900
      FLGFEQIFK                         40702000  48616000  37842000  40621000
      GANIASFVMVADAMLDQGDVF             56429000  86972000  50108000  54944000
      GCIISETGITSEQIHDIASAK              8370400   8934400   7776800   7904100
      GGLCVDLK                                NA        NA  13294000   9968700
      ICYAFMR                            9607500  12664000   8278600   8668400
      NSWEGVLTGK                        12479000  16368500  14892000  12989000
      SLEEIVDEYSTFSESK                   8899400  10091000   9655200   8573600
      VLPIVSVPER                        38299000  36292000  45293000  41889000
      VTISGSGNVAQYAALK                   4953700   3731485   1805400   3589830
      VTWENDNGEQEVAQGYR                  3413500   2769950   2269050   3205300
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK        NA        NA        NA        NA
      AANLGGVAVSGLEMAQNSQR              29672000  35306000  30967000  33135000
      ALVAQGVK                          12450000  14706000  19315000  12889000
      FIAEGSNMGSTPEAIAVFETAR            41313000  44409000  41112000  36661000
      GANIASFIK                         16999000  20105000  18688500  18302000
      GCIISETGITSEQVADISSAK              3749200   3172000   3574700   3610500
      HIGQDTDVPAGDIGVGGR                32530000  37002000  39590000  32416000
      IMINCFNECIDYAK                    14067000  13481000  12623000  13488000
      ITWTSER                           17868000  17810000  23153000  16760000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     88320000 102390000  93146000  90039000
      SLEQIVNEYSTFSENK                  56011000  50258000  59400000  53447000
      STATGPSEAVWYGPPK                  40218000  45023000  45462000  48656000
      VDIALPCATQNEVSGEEAK               49404000  51388000  54776000  54875000
      VIELGGTVVSLSDSK                   14980000  13746000  16469000  13965000
      VQYIAGARPWTHVQK                    5460050   3942750   3455200   7371750
      VTWENDKGEQEVAQGYR                 34850000  45930000  34888000  34159000
      AAGLTAAYAR                        60910000  63592000  75756000  55590000
      APEAEQVLSAAATFPIAQPATDVEAR        29710000  21845000  34877000  27494000
      AVQDNGESAFR                       12464000  15067000  15848000  13520000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     28071000  43997000  26263000  22185000
      GFTLAEVK                          29152000  34739000  30838000  30524000
      IAPRPLDLLRPVVR                    17729000  28612000  16477000  10253000
      IIVFPR                            48846000  44523000  57519000  55242000
      NQEIFDANVQR                       87065000  82847000  97118000  86931000
      TIGIAVDHR                         25122000  20284000  29831000  22671000
      VHFDQAGK                           9756600        NA  12641000  10128950
      VHFDQAGKK                               NA        NA        NA        NA
      ANELLINVK                        175640000 173790000 194190000 155330000
      ANGTTVLVGMPAGAK                         NA        NA        NA        NA
      ATDGGAHGVINVSVSEAAIEASTR         112140000  82689000 103380000  89717000
      CCSDVFNQVVK                       49834000  61137000  53120000  52898000
      DIVGAVLK                                NA  56266000  56161000  48921000
      EALDFFAR                         188800000 207980000 195260000 196550000
      EKDIVGAVLK                        39383000  35041000  36371000  30028000
      GVIFYESHGK                        94847000 120250000 106310000  86113000
      IGDYAGIK                         136670000 140680000 185230000 138400000
      LPLVGGHEGAGVVVGMGENVK            122650000 119390000 117670000 111850000
      SANLMAGHWVAISGAAGGLGSLAVQYAK      28729000  30301000  36144000  28932000
      SIGGEVFIDFTK                      66273000  77642000  67298000  69274000
      SIPETQK                           26807000  48695000  22957000  24072000
      SISIVGSYVGNR                     263390000 277410000 300530000 258680000
      VLGIDGGEGK                        21461000  11149000  32987000  19685000
      VLGIDGGEGKEELFR                  302640000 393230000 299540000 296920000
      VVGLSTLPEIYEK                    516890000 552670000 547410000 507720000
      YSGVCHTDLHAWHGDWPLPVK            100260000 107560000 106810000 146850000
      ANGTVVLVGLPAGAK                    4523400   5139700   5251000   5045100
      CSSDVFNHVVK                             NA   4415100        NA        NA
      DIPVPKPKPNELLINVK                       NA        NA        NA        NA
      VVGLSSLPEIYEK                     36993000  37893000  42723000  37048000
      DIPVPEPKPNEILINVK                 18878000  21317000  25243000  22329000
      EALDFFSR                           6149000   7301700   6106700   7275100
      GVIFYENK                           5329200   6083450   4909200   5064800
      IQQGTDLAEVAPILCAGVTVYK             8115600   8348700   7789300   7051200
      IVGLSELPK                         17116000  20323000  18659000  19479000
      NMVSDIQEATK                        5862500   4253700   6945900   6367200
      VLGIDAGEEK                              NA        NA        NA        NA
                                               5         6         7         8
      AAADALSDLEIK                            NA        NA        NA        NA
      AAADALSDLEIKDSK                   12612000   4766700        NA   3857100
      AEWALR                             5127000   5520400   5088000   5561300
      DEGLHTDFACLLFAHLK                       NA   4034700        NA        NA
      DIHDWNNR                           2031800   2201150        NA        NA
      ELETLREENR                         4106000   5888600        NA   4544500
      ESEFLFNAIHTIPEIGEK                26978000  39566000  32819000  33870000
      GMMPGLTFSNELICR                    7552500   8618700  10131000   8619700
      IVTEAVEIEQR                       17648000  27387000  16858000  18134000
      LLVAFGNK                           6784400  17607000  10856500   9171900
      LLVAFGNKK                               NA        NA        NA        NA
      NKPDPAIVEK                        18999000  25026000  16631000  20910000
      TNFFEK                             5853400   6233800   6488500   6357000
      TVLFPIK                           10083000  16620000  11783000  10174000
      VENPFDFMENISLAGK                  18359000  14936000  18244000  17782000
      WIQDADALFGER                      15309000  16929000  19946000  15877000
      YFLDALPVALLGMNADLMNQYVEFVADR      23777000  27010000  23475000  26459000
      AANLGGVAVSGLEMAQNSQK              10457000  10241000   9171700  10783000
      DAVWFGPPK                               NA        NA        NA        NA
      EIGYLFGAYR                        21641000        NA        NA        NA
      FHPSVNLSILK                        5435800   8597900   4928900   5049900
      FLGFEQIFK                         45495000  43813000  41632000  45829000
      GANIASFVMVADAMLDQGDVF             78097000  48322000  53524000  72863000
      GCIISETGITSEQIHDIASAK              8443900   8741000   7523500   8603100
      GGLCVDLK                          11138000  13324000  10058250  10731500
      ICYAFMR                           11127000  10013000   9772700  12064000
      NSWEGVLTGK                        14347000  16552000  12374000  14270000
      SLEEIVDEYSTFSESK                   8495100   6270500   8665300   8875100
      VLPIVSVPER                        34047000  47285000  44170000  42713000
      VTISGSGNVAQYAALK                   1694910   6089350   4744550        NA
      VTWENDNGEQEVAQGYR                  2857900        NA   3540200   2754200
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK        NA        NA        NA        NA
      AANLGGVAVSGLEMAQNSQR              35653000  27071000  30917000  35660000
      ALVAQGVK                          13377000  16599000  12316000  14473000
      FIAEGSNMGSTPEAIAVFETAR            41350000  33014000  39142000  42909000
      GANIASFIK                         18642000  21485000  18711000  22035000
      GCIISETGITSEQVADISSAK              3750800   4330000   4159300   3561200
      HIGQDTDVPAGDIGVGGR                32278000  48952000  31486000  36389000
      IMINCFNECIDYAK                    13703000  12631000  13346000  14127000
      ITWTSER                           16990000  19097000  18150000  20111000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR    100250000  11446000  81966000 100620000
      SLEQIVNEYSTFSENK                  54575000  60519000  58515000  52723000
      STATGPSEAVWYGPPK                  45139000  45577000  48340000  45932000
      VDIALPCATQNEVSGEEAK               47254000  56936000  53196000  56938000
      VIELGGTVVSLSDSK                   13320000  14777000  13868000  15646000
      VQYIAGARPWTHVQK                   10342000   4829150   9447700   7940200
      VTWENDKGEQEVAQGYR                 39544000  43698000  34957000  41604000
      AAGLTAAYAR                        59234000  73831000  56674000  68448000
      APEAEQVLSAAATFPIAQPATDVEAR        18843000  32425000  27790000  18780000
      AVQDNGESAFR                       11806000  15753000  14340000  14881000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     27642000  16185000  28015000  27163000
      GFTLAEVK                          31584000  35869000  28625000  35608000
      IAPRPLDLLRPVVR                    24360000  13232000  16928000  29313000
      IIVFPR                            44088000  62327000  55047000  50827000
      NQEIFDANVQR                       83524000  97828000  93185000  92385000
      TIGIAVDHR                         20934000  29048500  28133000  23862000
      VHFDQAGK                          10531000  24551000   9404800  11715000
      VHFDQAGKK                               NA        NA        NA        NA
      ANELLINVK                        159160000 239600000 165120000 170030000
      ANGTTVLVGMPAGAK                         NA  22228000  20300000        NA
      ATDGGAHGVINVSVSEAAIEASTR          75023000 138720000  89722000  96722000
      CCSDVFNQVVK                       56074000  56947000  49880000  53475000
      DIVGAVLK                          53999000  62136000        NA  51248000
      EALDFFAR                         186070000 215930000 198790000 215740000
      EKDIVGAVLK                        37416000  47384000  34168000  38143000
      GVIFYESHGK                        95943000 136440000  89159000 101610000
      IGDYAGIK                         130655000 169400000 139160000 131910000
      LPLVGGHEGAGVVVGMGENVK            114690000 159280000 128900000 149250000
      SANLMAGHWVAISGAAGGLGSLAVQYAK      29774000  22408000  22950000  37518000
      SIGGEVFIDFTK                      61694000  63157000  65128000  70890000
      SIPETQK                           25389000  11129000  33558000  21280000
      SISIVGSYVGNR                     255600000 336540000 272760000 288100000
      VLGIDGGEGK                        10732000  27786000  20564000  10019000
      VLGIDGGEGKEELFR                  331550000 366510000 301640000 364430000
      VVGLSTLPEIYEK                    514510000 563150000 542550000 555250000
      YSGVCHTDLHAWHGDWPLPVK            110030000 118290000 136520000 119170000
      ANGTVVLVGLPAGAK                    4766400   6666500   4201100   5176000
      CSSDVFNHVVK                        4047900   4399100   2461700   2850750
      DIPVPKPKPNELLINVK                       NA        NA        NA        NA
      VVGLSSLPEIYEK                     34603000  40097000  40341000  35832000
      DIPVPEPKPNEILINVK                 20059000  22053000  20546000  21474500
      EALDFFSR                           7094200   8147000   6789300   7007000
      GVIFYENK                           5928600   6391050        NA   5495200
      IQQGTDLAEVAPILCAGVTVYK             9297600   7746200   8002600   8705200
      IVGLSELPK                         18743000  22154000  18739000  22137000
      NMVSDIQEATK                        5760900   6084700   5901200   6515100
      VLGIDAGEEK                              NA        NA        NA        NA
                                               9
      AAADALSDLEIK                            NA
      AAADALSDLEIKDSK                         NA
      AEWALR                             5747200
      DEGLHTDFACLLFAHLK                  8548300
      DIHDWNNR                           2225200
      ELETLREENR                         8460600
      ESEFLFNAIHTIPEIGEK                36535000
      GMMPGLTFSNELICR                         NA
      IVTEAVEIEQR                       25508000
      LLVAFGNK                           7904400
      LLVAFGNKK                          2802500
      NKPDPAIVEK                        22252000
      TNFFEK                             6226000
      TVLFPIK                           10796000
      VENPFDFMENISLAGK                  15365000
      WIQDADALFGER                      15931000
      YFLDALPVALLGMNADLMNQYVEFVADR      31808000
      AANLGGVAVSGLEMAQNSQK              10396650
      DAVWFGPPK                               NA
      EIGYLFGAYR                        28783000
      FHPSVNLSILK                        8598000
      FLGFEQIFK                         61723000
      GANIASFVMVADAMLDQGDVF             99299000
      GCIISETGITSEQIHDIASAK             10779000
      GGLCVDLK                          12042000
      ICYAFMR                           10936000
      NSWEGVLTGK                        15028000
      SLEEIVDEYSTFSESK                   8548800
      VLPIVSVPER                        36401000
      VTISGSGNVAQYAALK                   3904350
      VTWENDNGEQEVAQGYR                  3463000
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK   6181500
      AANLGGVAVSGLEMAQNSQR              31639000
      ALVAQGVK                                NA
      FIAEGSNMGSTPEAIAVFETAR            43687000
      GANIASFIK                         18912000
      GCIISETGITSEQVADISSAK              6857850
      HIGQDTDVPAGDIGVGGR                39830000
      IMINCFNECIDYAK                    14142000
      ITWTSER                           19354000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR    103675000
      SLEQIVNEYSTFSENK                  48862000
      STATGPSEAVWYGPPK                  41778000
      VDIALPCATQNEVSGEEAK               48730000
      VIELGGTVVSLSDSK                   16526000
      VQYIAGARPWTHVQK                         NA
      VTWENDKGEQEVAQGYR                 54825000
      AAGLTAAYAR                        59012000
      APEAEQVLSAAATFPIAQPATDVEAR         3149500
      AVQDNGESAFR                       18130000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     32184000
      GFTLAEVK                          36671000
      IAPRPLDLLRPVVR                    10944000
      IIVFPR                             2242100
      NQEIFDANVQR                       84943000
      TIGIAVDHR                         32625000
      VHFDQAGK                          15879000
      VHFDQAGKK                          2097300
      ANELLINVK                        182820000
      ANGTTVLVGMPAGAK                   20851000
      ATDGGAHGVINVSVSEAAIEASTR         101040000
      CCSDVFNQVVK                       69073000
      DIVGAVLK                          60042000
      EALDFFAR                         212330000
      EKDIVGAVLK                        35459000
      GVIFYESHGK                       226950000
      IGDYAGIK                         129165000
      LPLVGGHEGAGVVVGMGENVK            128030000
      SANLMAGHWVAISGAAGGLGSLAVQYAK      46075000
      SIGGEVFIDFTK                      65134000
      SIPETQK                           33610000
      SISIVGSYVGNR                     271690000
      VLGIDGGEGK                         9843400
      VLGIDGGEGKEELFR                  460560000
      VVGLSTLPEIYEK                    557180000
      YSGVCHTDLHAWHGDWPLPVK            106120000
      ANGTVVLVGLPAGAK                    5874900
      CSSDVFNHVVK                        5483200
      DIPVPKPKPNELLINVK                  6594800
      VVGLSSLPEIYEK                     34226000
      DIPVPEPKPNEILINVK                 21048000
      EALDFFSR                           6971000
      GVIFYENK                           7300500
      IQQGTDLAEVAPILCAGVTVYK             8480300
      IVGLSELPK                         17048000
      NMVSDIQEATK                        7612200
      VLGIDAGEEK                              NA

---

    Code
      as.data.frame(tail(SummarizedExperiment::rowData(D2), n = 1000))
    Output
                                                               Sequence
      AAADALSDLEIK                                         AAADALSDLEIK
      AAADALSDLEIKDSK                                   AAADALSDLEIKDSK
      AEWALR                                                     AEWALR
      DEGLHTDFACLLFAHLK                               DEGLHTDFACLLFAHLK
      DIHDWNNR                                                 DIHDWNNR
      ELETLREENR                                             ELETLREENR
      ESEFLFNAIHTIPEIGEK                             ESEFLFNAIHTIPEIGEK
      GMMPGLTFSNELICR                                   GMMPGLTFSNELICR
      IVTEAVEIEQR                                           IVTEAVEIEQR
      LLVAFGNK                                                 LLVAFGNK
      LLVAFGNKK                                               LLVAFGNKK
      NKPDPAIVEK                                             NKPDPAIVEK
      TNFFEK                                                     TNFFEK
      TVLFPIK                                                   TVLFPIK
      VENPFDFMENISLAGK                                 VENPFDFMENISLAGK
      WIQDADALFGER                                         WIQDADALFGER
      YFLDALPVALLGMNADLMNQYVEFVADR         YFLDALPVALLGMNADLMNQYVEFVADR
      AANLGGVAVSGLEMAQNSQK                         AANLGGVAVSGLEMAQNSQK
      DAVWFGPPK                                               DAVWFGPPK
      EIGYLFGAYR                                             EIGYLFGAYR
      FHPSVNLSILK                                           FHPSVNLSILK
      FLGFEQIFK                                               FLGFEQIFK
      GANIASFVMVADAMLDQGDVF                       GANIASFVMVADAMLDQGDVF
      GCIISETGITSEQIHDIASAK                       GCIISETGITSEQIHDIASAK
      GGLCVDLK                                                 GGLCVDLK
      ICYAFMR                                                   ICYAFMR
      NSWEGVLTGK                                             NSWEGVLTGK
      SLEEIVDEYSTFSESK                                 SLEEIVDEYSTFSESK
      VLPIVSVPER                                             VLPIVSVPER
      VTISGSGNVAQYAALK                                 VTISGSGNVAQYAALK
      VTWENDNGEQEVAQGYR                               VTWENDNGEQEVAQGYR
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK YVAGARPWTHVSNVDIALPCATQNEVSGDEAK
      AANLGGVAVSGLEMAQNSQR                         AANLGGVAVSGLEMAQNSQR
      ALVAQGVK                                                 ALVAQGVK
      FIAEGSNMGSTPEAIAVFETAR                     FIAEGSNMGSTPEAIAVFETAR
      GANIASFIK                                               GANIASFIK
      GCIISETGITSEQVADISSAK                       GCIISETGITSEQVADISSAK
      HIGQDTDVPAGDIGVGGR                             HIGQDTDVPAGDIGVGGR
      IMINCFNECIDYAK                                     IMINCFNECIDYAK
      ITWTSER                                                   ITWTSER
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR       SEPEFQQAYEEVVSSLEDSTLFEQHPEYR
      SLEQIVNEYSTFSENK                                 SLEQIVNEYSTFSENK
      STATGPSEAVWYGPPK                                 STATGPSEAVWYGPPK
      VDIALPCATQNEVSGEEAK                           VDIALPCATQNEVSGEEAK
      VIELGGTVVSLSDSK                                   VIELGGTVVSLSDSK
      VQYIAGARPWTHVQK                                   VQYIAGARPWTHVQK
      VTWENDKGEQEVAQGYR                               VTWENDKGEQEVAQGYR
      AAGLTAAYAR                                             AAGLTAAYAR
      APEAEQVLSAAATFPIAQPATDVEAR             APEAEQVLSAAATFPIAQPATDVEAR
      AVQDNGESAFR                                           AVQDNGESAFR
      DGKAPEAEQVLSAAATFPIAQPATDVEAR       DGKAPEAEQVLSAAATFPIAQPATDVEAR
      GFTLAEVK                                                 GFTLAEVK
      IAPRPLDLLRPVVR                                     IAPRPLDLLRPVVR
      IIVFPR                                                     IIVFPR
      NQEIFDANVQR                                           NQEIFDANVQR
      TIGIAVDHR                                               TIGIAVDHR
      VHFDQAGK                                                 VHFDQAGK
      VHFDQAGKK                                               VHFDQAGKK
      ANELLINVK                                               ANELLINVK
      ANGTTVLVGMPAGAK                                   ANGTTVLVGMPAGAK
      ATDGGAHGVINVSVSEAAIEASTR                 ATDGGAHGVINVSVSEAAIEASTR
      CCSDVFNQVVK                                           CCSDVFNQVVK
      DIVGAVLK                                                 DIVGAVLK
      EALDFFAR                                                 EALDFFAR
      EKDIVGAVLK                                             EKDIVGAVLK
      GVIFYESHGK                                             GVIFYESHGK
      IGDYAGIK                                                 IGDYAGIK
      LPLVGGHEGAGVVVGMGENVK                       LPLVGGHEGAGVVVGMGENVK
      SANLMAGHWVAISGAAGGLGSLAVQYAK         SANLMAGHWVAISGAAGGLGSLAVQYAK
      SIGGEVFIDFTK                                         SIGGEVFIDFTK
      SIPETQK                                                   SIPETQK
      SISIVGSYVGNR                                         SISIVGSYVGNR
      VLGIDGGEGK                                             VLGIDGGEGK
      VLGIDGGEGKEELFR                                   VLGIDGGEGKEELFR
      VVGLSTLPEIYEK                                       VVGLSTLPEIYEK
      YSGVCHTDLHAWHGDWPLPVK                       YSGVCHTDLHAWHGDWPLPVK
      ANGTVVLVGLPAGAK                                   ANGTVVLVGLPAGAK
      CSSDVFNHVVK                                           CSSDVFNHVVK
      DIPVPKPKPNELLINVK                               DIPVPKPKPNELLINVK
      VVGLSSLPEIYEK                                       VVGLSSLPEIYEK
      DIPVPEPKPNEILINVK                               DIPVPEPKPNEILINVK
      EALDFFSR                                                 EALDFFSR
      GVIFYENK                                                 GVIFYENK
      IQQGTDLAEVAPILCAGVTVYK                     IQQGTDLAEVAPILCAGVTVYK
      IVGLSELPK                                               IVGLSELPK
      NMVSDIQEATK                                           NMVSDIQEATK
      VLGIDAGEEK                                             VLGIDAGEEK

---

    Code
      as.data.frame(tail(SummarizedExperiment::colData(D2), n = 1000))
    Output
        group
      1     1
      2     2
      3     3
      4     4
      5     5
      6     6
      7     7
      8     8
      9     9

# test calculatePeptideRatios

    Code
      D1
    Output
      class: SummarizedExperiment 
      dim: 87 36 
      metadata(0):
      assays(1): logRatios
      rownames(87): AAADALSDLEIK AAADALSDLEIKDSK ... NMVSDIQEATK VLGIDAGEEK
      rowData names(1): Sequence
      colnames(36): logRatio_1_2 logRatio_1_3 ... logRatio_7_9 logRatio_8_9
      colData names(1): comparison

---

    Code
      SummarizedExperiment::assays(D1)$logRatios
    Output
                                       logRatio_1_2 logRatio_1_3  logRatio_1_4
      AAADALSDLEIK                               NA           NA            NA
      AAADALSDLEIKDSK                   0.990798156  0.206301207  0.0638614470
      AEWALR                            0.277627128  0.173103910  0.0201120494
      DEGLHTDFACLLFAHLK                 0.520190758  0.774659104  0.6540415735
      DIHDWNNR                                   NA           NA            NA
      ELETLREENR                                 NA           NA            NA
      ESEFLFNAIHTIPEIGEK                0.095884188  0.008068653 -0.0992333027
      GMMPGLTFSNELICR                  -0.086926627  0.063867600  0.1489593812
      IVTEAVEIEQR                       0.005326543  0.255680238 -0.0872181736
      LLVAFGNK                         -0.689461662  0.163538124 -0.2434747178
      LLVAFGNKK                                  NA           NA            NA
      NKPDPAIVEK                        0.184692807  0.138061650 -0.1086336718
      TNFFEK                           -0.212800567  0.132000539 -0.0625490368
      TVLFPIK                          -0.254663486  0.095447879 -0.2573261825
      VENPFDFMENISLAGK                 -0.032062967 -0.057916143 -0.1775071422
      WIQDADALFGER                     -0.081982196  0.075858794  0.1081800806
      YFLDALPVALLGMNADLMNQYVEFVADR     -0.239158120 -0.073719602 -0.2245808240
      AANLGGVAVSGLEMAQNSQK              0.411957645  0.388555899  0.3660073434
      DAVWFGPPK                                  NA           NA            NA
      EIGYLFGAYR                                 NA           NA            NA
      FHPSVNLSILK                                NA           NA            NA
      FLGFEQIFK                         0.291243201 -0.053798388 -0.0087176454
      GANIASFVMVADAMLDQGDVF             0.592075811 -0.173684216 -0.0504132540
      GCIISETGITSEQIHDIASAK             0.031837183 -0.027777543 -0.1562673771
      GGLCVDLK                                   NA           NA            NA
      ICYAFMR                           0.427870687 -0.148057613 -0.0910838438
      NSWEGVLTGK                                 NA  0.224504280 -0.0359494109
      SLEEIVDEYSTFSESK                  0.103271399  0.006480548  0.0005581076
      VLPIVSVPER                       -0.083933176  0.287482565  0.0873369792
      VTISGSGNVAQYAALK                           NA -0.192193173            NA
      VTWENDNGEQEVAQGYR                          NA           NA            NA
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA            NA
      AANLGGVAVSGLEMAQNSQR             -0.114393489  0.131127888  0.1500252445
      ALVAQGVK                          0.263318925  0.465405316            NA
      FIAEGSNMGSTPEAIAVFETAR            0.297536245  0.243090129  0.0234306072
      GANIASFIK                         0.188506076           NA  0.0421217281
      GCIISETGITSEQVADISSAK            -0.243257398 -0.099765989 -0.0511618030
      HIGQDTDVPAGDIGVGGR                0.206427446  0.192400561  0.0103538415
      IMINCFNECIDYAK                   -0.043397418 -0.139836016 -0.0224264962
      ITWTSER                          -0.568146142  0.244476782            NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     0.266623546  0.495226033  0.4931219773
      SLEQIVNEYSTFSENK                 -0.267719116 -0.048247619 -0.1371977634
      STATGPSEAVWYGPPK                  0.056385104  0.153980712  0.1483339925
      VDIALPCATQNEVSGEEAK               0.118489937  0.192708642  0.1856070549
      VIELGGTVVSLSDSK                  -0.116937083  0.098869710 -0.1078708986
      VQYIAGARPWTHVQK                            NA           NA            NA
      VTWENDKGEQEVAQGYR                 0.367863441 -0.030884983 -0.0684876430
      AAGLTAAYAR                        0.096709152  0.239896021 -0.0963121297
      APEAEQVLSAAATFPIAQPATDVEAR       -0.351620711  0.207397368 -0.0713328801
      AVQDNGESAFR                       0.288672226  0.236872833  0.1217255528
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     0.630849081 -0.230092936 -0.2336273152
      GFTLAEVK                          0.233406229  0.083465068  0.0699820572
      IAPRPLDLLRPVVR                    0.626732223  0.002673044 -0.4841394054
      IIVFPR                           -0.668130351  0.241663598  0.1688064074
      NQEIFDANVQR                      -0.045836459  0.110545839  0.0055819374
      TIGIAVDHR                        -0.337004859  0.178992085 -0.1371829859
      VHFDQAGK                                   NA  0.272597928            NA
      VHFDQAGKK                                  NA           NA            NA
      ANELLINVK                         0.017405583  0.172551888 -0.1940217980
      ANGTTVLVGMPAGAK                            NA           NA            NA
      ATDGGAHGVINVSVSEAAIEASTR         -0.536764290 -0.159707871 -0.4349050003
      CCSDVFNQVVK                       0.269715364  0.083133758  0.0841351674
      DIVGAVLK                                   NA           NA            NA
      EALDFFAR                          0.132105131  0.144830887  0.0542332530
      EKDIVGAVLK                       -0.169943302 -0.121676214 -0.2901683119
      GVIFYESHGK                        0.316442537  0.095436486 -0.1186754364
      IGDYAGIK                          0.039652046  0.353006221 -0.0094654244
      LPLVGGHEGAGVVVGMGENVK            -0.114630446 -0.114178044 -0.0672918909
      SANLMAGHWVAISGAAGGLGSLAVQYAK      0.107345376  0.477540651  0.0972193216
      SIGGEVFIDFTK                      0.221763392  0.104398497  0.0050654952
      SIPETQK                           0.603409892 -0.388440003  0.0384623311
      SISIVGSYVGNR                      0.047614392  0.157469806 -0.0417438465
      VLGIDGGEGK                                 NA           NA            NA
      VLGIDGGEGKEELFR                   0.336461526 -0.013954428  0.0029772749
      VVGLSTLPEIYEK                     0.082929261  0.087019813 -0.0169171357
      YSGVCHTDLHAWHGDWPLPVK             0.128965469  0.148478628  0.5734287263
      ANGTVVLVGLPAGAK                   0.227487450  0.121933395 -0.0002154938
      CSSDVFNHVVK                                NA           NA            NA
      DIPVPKPKPNELLINVK                          NA           NA            NA
      VVGLSSLPEIYEK                    -0.005253044  0.189987639 -0.0577421214
      DIPVPEPKPNEILINVK                 0.145259189  0.359725097  0.2392425082
      EALDFFSR                          0.135391496  0.011997867  0.1930549930
      GVIFYENK                                   NA  0.015226332 -0.0107930991
      IQQGTDLAEVAPILCAGVTVYK            0.151852213  0.048403367 -0.0128978133
      IVGLSELPK                         0.245356505  0.135029250  0.1555402567
      NMVSDIQEATK                                NA  0.139035844  0.0735710028
      VLGIDAGEEK                                 NA           NA            NA
                                        logRatio_1_5 logRatio_1_6  logRatio_1_7
      AAADALSDLEIK                                NA           NA            NA
      AAADALSDLEIKDSK                   0.8224705863 -0.115838794            NA
      AEWALR                            0.1332858376  0.146864120  0.0888834668
      DEGLHTDFACLLFAHLK                           NA           NA            NA
      DIHDWNNR                                    NA           NA            NA
      ELETLREENR                                  NA           NA            NA
      ESEFLFNAIHTIPEIGEK               -0.2130617746  0.168076244 -0.1053327437
      GMMPGLTFSNELICR                             NA  0.038616681  0.3271444822
      IVTEAVEIEQR                      -0.0803997817  0.481377893 -0.1699717579
      LLVAFGNK                         -0.8217584022  0.477057105            NA
      LLVAFGNKK                                   NA           NA            NA
      NKPDPAIVEK                        0.1175762401  0.557508973 -0.0225588343
      TNFFEK                           -0.1204108509 -0.040512503  0.0479617672
      TVLFPIK                          -0.3990951775  0.296318170 -0.2336618450
      VENPFDFMENISLAGK                 -0.0382084716 -0.304851514 -0.0292927021
      WIQDADALFGER                     -0.2566439480 -0.165715197  0.1211555029
      YFLDALPVALLGMNADLMNQYVEFVADR     -0.2583294989 -0.054509326 -0.2908246924
      AANLGGVAVSGLEMAQNSQK              0.2354930187  0.426771689  0.0992330969
      DAVWFGPPK                                   NA           NA            NA
      EIGYLFGAYR                                  NA           NA            NA
      FHPSVNLSILK                                 NA           NA            NA
      FLGFEQIFK                         0.1553556145  0.102213611  0.0155672225
      GANIASFVMVADAMLDQGDVF             0.4772334161 -0.237455528 -0.1000384394
      GCIISETGITSEQIHDIASAK            -0.0395244111  0.124461720 -0.2040834286
      GGLCVDLK                                    NA           NA            NA
      ICYAFMR                           0.2135290916  0.061449437 -0.0098807222
      NSWEGVLTGK                        0.0593788453  0.354199654 -0.1254636389
      SLEEIVDEYSTFSESK                 -0.1304741184 -0.434444137 -0.0852141252
      VLPIVSVPER                       -0.1317628796  0.255484999  0.1789155173
      VTISGSGNVAQYAALK                            NA           NA            NA
      VTWENDNGEQEVAQGYR                           NA           NA            NA
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK            NA           NA            NA
      AANLGGVAVSGLEMAQNSQR              0.3139613770 -0.022870408 -0.2265358584
      ALVAQGVK                          0.0738938722  0.437884586  0.0101304581
      FIAEGSNMGSTPEAIAVFETAR            0.0562452862 -0.242343453  0.1215673630
      GANIASFIK                         0.0696305993  0.307395569  0.0538434743
      GCIISETGITSEQVADISSAK             0.0562628041  0.337005256  0.0456504841
      HIGQDTDVPAGDIGVGGR                0.1046092046  0.564313413 -0.0330933038
      IMINCFNECIDYAK                   -0.0001404049  0.036460063 -0.0290009407
      ITWTSER                          -0.1072550891           NA  0.0561692197
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     0.6542949357 -0.819666142  0.4127430495
      SLEQIVNEYSTFSENK                 -0.1115848795 -0.043459030 -0.0333245130
      STATGPSEAVWYGPPK                 -0.0016575095           NA  0.1367007948
      VDIALPCATQNEVSGEEAK              -0.0254804303  0.217498302  0.1207751583
      VIELGGTVVSLSDSK                  -0.1641108952  0.009094154 -0.1174578674
      VQYIAGARPWTHVQK                             NA           NA            NA
      VTWENDKGEQEVAQGYR                 0.1637676167  0.265610501  0.0053389100
      AAGLTAAYAR                        0.0053535651  0.320580502 -0.0384416912
      APEAEQVLSAAATFPIAQPATDVEAR       -0.7164521598  0.136951159 -0.0658919252
      AVQDNGESAFR                      -0.4846211632  0.353023068  0.1584469296
      DGKAPEAEQVLSAAATFPIAQPATDVEAR    -0.0215680643 -0.890066158  0.0478288947
      GFTLAEVK                          0.1391517631  0.265465436  0.0009553399
      IAPRPLDLLRPVVR                    0.1997644381 -0.157089596 -0.2319524524
      IIVFPR                           -0.1391772150  0.309666338  0.1430743603
      NQEIFDANVQR                      -0.0506729193  0.165431656  0.1124625893
      TIGIAVDHR                        -0.2185271882           NA  0.1099800195
      VHFDQAGK                          0.0458441836  1.209499019 -0.1412524046
      VHFDQAGKK                                   NA           NA            NA
      ANELLINVK                        -0.0822853088  0.489818766 -0.0678567718
      ANGTTVLVGMPAGAK                             NA           NA            NA
      ATDGGAHGVINVSVSEAAIEASTR         -0.6388642744  0.312932345 -0.3383822471
      CCSDVFNQVVK                       0.1458712736  0.221365649  0.0158593776
      DIVGAVLK                                    NA           NA            NA
      EALDFFAR                          0.0236363324  0.201246647  0.0685376946
      EKDIVGAVLK                       -0.1324017291  0.246514122 -0.2315085665
      GVIFYESHGK                        0.0201700105  0.543766679 -0.0787634302
      IGDYAGIK                                    NA  0.312051567  0.0223558981
      LPLVGGHEGAGVVVGMGENVK            -0.2535462160  0.229411607 -0.0069945814
      SANLMAGHWVAISGAAGGLGSLAVQYAK      0.0813563913  0.013510918 -0.3829049557
      SIGGEVFIDFTK                     -0.1401273114 -0.090106793 -0.0342322608
      SIPETQK                          -0.1236809801 -1.193022612  0.0552311424
      SISIVGSYVGNR                     -0.0896520004  0.328498586  0.0334203876
      VLGIDGGEGK                                  NA           NA            NA
      VLGIDGGEGKEELFR                   0.1120085781  0.255325178 -0.0043343518
      VVGLSTLPEIYEK                    -0.0163369245  0.094169990  0.0568592533
      YSGVCHTDLHAWHGDWPLPVK             0.1165588109  0.219819151  0.4708742099
      ANGTVVLVGLPAGAK                   0.0732975424  0.336811369 -0.0840576979
      CSSDVFNHVVK                                 NA           NA            NA
      DIPVPKPKPNELLINVK                           NA           NA            NA
      VVGLSSLPEIYEK                    -0.1220624212  0.077567822  0.0718794614
      DIPVPEPKPNEILINVK                -0.4467566062  0.218185906  0.1056013107
      EALDFFSR                                    NA           NA  0.0336353764
      GVIFYENK                                    NA           NA            NA
      IQQGTDLAEVAPILCAGVTVYK            0.1798689432  0.010753986  0.1436694400
      IVGLSELPK                         0.0996951424  0.376210689  0.1318959252
      NMVSDIQEATK                      -0.0289518543 -0.049872968  0.0147733741
      VLGIDAGEEK                                  NA           NA            NA
                                       logRatio_1_8 logRatio_1_9  logRatio_2_3
      AAADALSDLEIK                               NA           NA            NA
      AAADALSDLEIKDSK                  -0.007158832           NA -0.7844969486
      AEWALR                            0.159892466  0.231620482 -0.1045232183
      DEGLHTDFACLLFAHLK                          NA  1.712775676  0.2544683453
      DIHDWNNR                                   NA           NA            NA
      ELETLREENR                                 NA           NA            NA
      ESEFLFNAIHTIPEIGEK                0.001155431  0.014797286 -0.0878155343
      GMMPGLTFSNELICR                   0.172672239           NA  0.1507942270
      IVTEAVEIEQR                      -0.042359343  0.453805440  0.2503536948
      LLVAFGNK                         -0.468867263 -0.659305212  0.8529997858
      LLVAFGNKK                                  NA           NA            NA
      NKPDPAIVEK                        0.279713364  0.309246142 -0.0466311575
      TNFFEK                           -0.006530520  0.021005839  0.3448011061
      TVLFPIK                          -0.303709964 -0.298390061  0.3501113653
      VENPFDFMENISLAGK                 -0.091502235 -0.268149286 -0.0258531758
      WIQDADALFGER                     -0.177145021 -0.219425584  0.1578409900
      YFLDALPVALLGMNADLMNQYVEFVADR     -0.106183380  0.144000871  0.1654385188
      AANLGGVAVSGLEMAQNSQK              0.471090313           NA -0.0234017463
      DAVWFGPPK                                  NA           NA            NA
      EIGYLFGAYR                                 NA           NA            NA
      FHPSVNLSILK                                NA           NA  0.4483464726
      FLGFEQIFK                         0.193126379  0.601279045 -0.3450415895
      GANIASFVMVADAMLDQGDVF             0.358079891  0.933342124 -0.7657600275
      GCIISETGITSEQIHDIASAK            -0.013434029  0.272120371 -0.0596147258
      GGLCVDLK                                   NA           NA            NA
      ICYAFMR                           0.369583144  0.176452109 -0.5759282996
      NSWEGVLTGK                        0.125016088  0.158869560            NA
      SLEEIVDEYSTFSESK                  0.017583152 -0.080125931 -0.0967908508
      VLPIVSVPER                        0.167332389 -0.126717517  0.3714157412
      VTISGSGNVAQYAALK                           NA           NA            NA
      VTWENDNGEQEVAQGYR                          NA           NA            NA
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA            NA
      AANLGGVAVSGLEMAQNSQR             -0.102120867  0.083725758  0.2455213770
      ALVAQGVK                          0.205268227           NA  0.2020863903
      FIAEGSNMGSTPEAIAVFETAR            0.244305390  0.271580258 -0.0544461163
      GANIASFIK                         0.338062748  0.086538441            NA
      GCIISETGITSEQVADISSAK             0.053769806           NA  0.1434914085
      HIGQDTDVPAGDIGVGGR                0.277022353  0.399095701 -0.0140268853
      IMINCFNECIDYAK                    0.034850657  0.082932027 -0.0964385972
      ITWTSER                           0.109482298           NA  0.8126229243
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     0.665863618           NA  0.2286024875
      SLEQIVNEYSTFSENK                 -0.152371298 -0.675817348  0.2194714967
      STATGPSEAVWYGPPK                  0.128445442 -0.003183211  0.0975956083
      VDIALPCATQNEVSGEEAK               0.258193646 -0.426913270  0.0742187056
      VIELGGTVVSLSDSK                   0.017245981  0.208174125  0.2158067925
      VQYIAGARPWTHVQK                            NA           NA            NA
      VTWENDKGEQEVAQGYR                 0.235906652  0.716837780 -0.3987484248
      AAGLTAAYAR                        0.203518660 -0.009590376  0.1431868684
      APEAEQVLSAAATFPIAQPATDVEAR       -0.712303540 -2.936524553  0.5590180790
      AVQDNGESAFR                       0.246292057  0.501448147 -0.0517993927
      DGKAPEAEQVLSAAATFPIAQPATDVEAR    -0.051776102  0.099132063 -0.8609420166
      GFTLAEVK                          0.239880897  0.316556546 -0.1499411609
      IAPRPLDLLRPVVR                    0.669464542 -0.335219591 -0.6240591786
      IIVFPR                            0.046846039 -1.553596893  0.9097939495
      NQEIFDANVQR                       0.122967122 -0.035490935  0.1563822987
      TIGIAVDHR                        -0.063168764  0.336214238  0.5159969434
      VHFDQAGK                          0.156692835  0.556758911            NA
      VHFDQAGKK                                  NA           NA            NA
      ANELLINVK                        -0.036795730  0.076999944  0.1551463049
      ANGTTVLVGMPAGAK                            NA           NA            NA
      ATDGGAHGVINVSVSEAAIEASTR         -0.264332021 -0.112006950  0.3770564188
      CCSDVFNQVVK                       0.119183951  0.485627088 -0.1865816056
      DIVGAVLK                                   NA           NA            NA
      EALDFFAR                          0.186008693  0.200691717  0.0127257564
      EKDIVGAVLK                       -0.084590938  0.156611208  0.0482670884
      GVIFYESHGK                        0.086175784  1.137496163 -0.2210060502
      IGDYAGIK                         -0.049745006           NA  0.3133541758
      LPLVGGHEGAGVVVGMGENVK             0.133549648 -0.046025181  0.0004524024
      SANLMAGHWVAISGAAGGLGSLAVQYAK      0.412849540  0.890315663  0.3701952750
      SIGGEVFIDFTK                      0.100476875  0.027711505 -0.1173648949
      SIPETQK                          -0.321194283  0.584887660 -0.9918498954
      SISIVGSYVGNR                      0.098159288  0.063896024  0.1098554138
      VLGIDGGEGK                                 NA           NA  1.3346718183
      VLGIDGGEGKEELFR                   0.238116313  0.623941907 -0.3504159543
      VVGLSTLPEIYEK                     0.133739343  0.079264966  0.0040905525
      YSGVCHTDLHAWHGDWPLPVK             0.209923833  0.095251996  0.0195131594
      ANGTVVLVGLPAGAK                   0.193237520  0.338654436 -0.1055540551
      CSSDVFNHVVK                                NA           NA            NA
      DIPVPKPKPNELLINVK                          NA           NA            NA
      VVGLSSLPEIYEK                    -0.088344942 -0.158903942  0.1952406833
      DIPVPEPKPNEILINVK                          NA -0.105196990  0.2144659083
      EALDFFSR                          0.301009593  0.090298305 -0.1233936285
      GVIFYENK                          0.080037978  0.440906884            NA
      IQQGTDLAEVAPILCAGVTVYK            0.223830125  0.155379843 -0.1034488468
      IVGLSELPK                         0.313893387 -0.034532726 -0.1103272553
      NMVSDIQEATK                       0.004836273  0.326806899            NA
      VLGIDAGEEK                                 NA           NA            NA
                                       logRatio_2_4 logRatio_2_5 logRatio_2_6
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                  -0.926936709 -0.168327570  -1.10663695
      AEWALR                           -0.257515079 -0.144341291  -0.13076301
      DEGLHTDFACLLFAHLK                 0.133850815           NA           NA
      DIHDWNNR                                   NA           NA           NA
      ELETLREENR                                 NA -0.514541378   0.01300647
      ESEFLFNAIHTIPEIGEK               -0.195117490 -0.308945962   0.07219206
      GMMPGLTFSNELICR                   0.235886008           NA   0.12554331
      IVTEAVEIEQR                      -0.092544717 -0.085726325   0.47605135
      LLVAFGNK                          0.445986944 -0.132296740   1.16651877
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                       -0.293326479 -0.067116567   0.37281617
      TNFFEK                            0.150251530  0.092389716   0.17228806
      TVLFPIK                          -0.002662696 -0.144431692   0.55098166
      VENPFDFMENISLAGK                 -0.145444175 -0.006145504  -0.27278855
      WIQDADALFGER                      0.190162277 -0.174661752  -0.08373300
      YFLDALPVALLGMNADLMNQYVEFVADR      0.014577296 -0.019171378   0.18464879
      AANLGGVAVSGLEMAQNSQK             -0.045950301 -0.176464626   0.01481404
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                       0.211115535  0.114155217   0.92229585
      FLGFEQIFK                        -0.299960847 -0.135887587  -0.18902959
      GANIASFVMVADAMLDQGDVF            -0.642489065 -0.114842395  -0.82953134
      GCIISETGITSEQIHDIASAK            -0.188104560 -0.071361594   0.09262454
      GGLCVDLK                                   NA           NA           NA
      ICYAFMR                          -0.518954531 -0.214341595  -0.36642125
      NSWEGVLTGK                                 NA           NA           NA
      SLEEIVDEYSTFSESK                 -0.102713291 -0.233745517  -0.53771554
      VLPIVSVPER                        0.171270155 -0.047829703   0.33941817
      VTISGSGNVAQYAALK                           NA           NA           NA
      VTWENDNGEQEVAQGYR                          NA           NA           NA
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR              0.264418733  0.428354866   0.09152308
      ALVAQGVK                                   NA -0.189425053   0.17456566
      FIAEGSNMGSTPEAIAVFETAR           -0.274105638 -0.241290959  -0.53987970
      GANIASFIK                        -0.146384348 -0.118875477   0.11888949
      GCIISETGITSEQVADISSAK             0.192095595  0.299520202   0.58026265
      HIGQDTDVPAGDIGVGGR               -0.196073605 -0.101818242   0.35788597
      IMINCFNECIDYAK                    0.020970922  0.043257013   0.07985748
      ITWTSER                                    NA  0.460891053           NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     0.226498432  0.387671390  -1.08628969
      SLEQIVNEYSTFSENK                  0.130521352  0.156134236   0.22426009
      STATGPSEAVWYGPPK                  0.091948888 -0.058042614           NA
      VDIALPCATQNEVSGEEAK               0.067117118 -0.143970367   0.09900837
      VIELGGTVVSLSDSK                   0.009066184 -0.047173813   0.12603124
      VQYIAGARPWTHVQK                            NA           NA           NA
      VTWENDKGEQEVAQGYR                -0.436351084 -0.204095825  -0.10225294
      AAGLTAAYAR                       -0.193021282 -0.091355587   0.22387135
      APEAEQVLSAAATFPIAQPATDVEAR        0.280287831 -0.364831449   0.48857187
      AVQDNGESAFR                      -0.166946673 -0.773293389   0.06435084
      DGKAPEAEQVLSAAATFPIAQPATDVEAR    -0.864476396 -0.652417145  -1.52091524
      GFTLAEVK                         -0.163424171 -0.094254465   0.03205921
      IAPRPLDLLRPVVR                   -1.110871628 -0.426967785  -0.78382182
      IIVFPR                            0.836936759  0.528953136   0.97779669
      NQEIFDANVQR                       0.051418397 -0.004836460   0.21126812
      TIGIAVDHR                         0.199821873  0.118477670           NA
      VHFDQAGK                                   NA           NA           NA
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                        -0.211427381 -0.099690892   0.47241318
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR          0.101859290 -0.102099984   0.84969664
      CCSDVFNQVVK                      -0.185580197 -0.123844090  -0.04834972
      DIVGAVLK                                   NA           NA           NA
      EALDFFAR                         -0.077871878 -0.108468798   0.06914152
      EKDIVGAVLK                       -0.120225009  0.037541573   0.41645742
      GVIFYESHGK                       -0.435117973 -0.296272526   0.22732414
      IGDYAGIK                         -0.049117470           NA   0.27239952
      LPLVGGHEGAGVVVGMGENVK             0.047338555 -0.138915770   0.34404205
      SANLMAGHWVAISGAAGGLGSLAVQYAK     -0.010126054 -0.025988984  -0.09383446
      SIGGEVFIDFTK                     -0.216697897 -0.361890704  -0.31187019
      SIPETQK                          -0.564947561 -0.727090872  -1.79643250
      SISIVGSYVGNR                     -0.089358239 -0.137266393   0.28088419
      VLGIDGGEGK                        0.812120173 -0.113282325   1.30903679
      VLGIDGGEGKEELFR                  -0.333484251 -0.224452948  -0.08113635
      VVGLSTLPEIYEK                    -0.099846396 -0.099266185   0.01124073
      YSGVCHTDLHAWHGDWPLPVK             0.444463257 -0.012406658   0.09085368
      ANGTVVLVGLPAGAK                  -0.227702944 -0.154189908   0.10932392
      CSSDVFNHVVK                                NA           NA  -0.06025549
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                    -0.052489077 -0.116809377   0.08282087
      DIPVPEPKPNEILINVK                 0.093983319 -0.592015795   0.07292672
      EALDFFSR                          0.057663497           NA           NA
      GVIFYENK                                   NA           NA           NA
      IQQGTDLAEVAPILCAGVTVYK           -0.164750027  0.028016730  -0.14109823
      IVGLSELPK                        -0.089816249 -0.145661363   0.13085418
      NMVSDIQEATK                                NA           NA           NA
      VLGIDAGEEK                                 NA           NA           NA
                                        logRatio_2_7 logRatio_2_8 logRatio_2_9
      AAADALSDLEIK                                NA           NA           NA
      AAADALSDLEIKDSK                             NA -0.997956988           NA
      AEWALR                           -0.1887436613 -0.117734662 -0.046006646
      DEGLHTDFACLLFAHLK                           NA           NA  1.192584918
      DIHDWNNR                                    NA           NA           NA
      ELETLREENR                                  NA -0.299306037  0.506448217
      ESEFLFNAIHTIPEIGEK               -0.2012169312 -0.094728756 -0.081086901
      GMMPGLTFSNELICR                   0.4140711088  0.259598866           NA
      IVTEAVEIEQR                      -0.1752983010 -0.047685886  0.448478897
      LLVAFGNK                                    NA  0.220594399  0.030156450
      LLVAFGNKK                                   NA           NA  0.590085766
      NKPDPAIVEK                       -0.2072516414  0.095020557  0.124553335
      TNFFEK                            0.2607623338  0.206270047  0.233806405
      TVLFPIK                           0.0210016411 -0.049046478 -0.043726574
      VENPFDFMENISLAGK                  0.0027702651 -0.059439268 -0.236086319
      WIQDADALFGER                      0.2031376994 -0.095162824 -0.137443388
      YFLDALPVALLGMNADLMNQYVEFVADR     -0.0516665720  0.132974740  0.383158991
      AANLGGVAVSGLEMAQNSQK             -0.3127245480  0.059132668           NA
      DAVWFGPPK                                   NA           NA           NA
      EIGYLFGAYR                                  NA           NA  0.155593414
      FHPSVNLSILK                      -0.0230942179 -0.020527604  1.007246594
      FLGFEQIFK                        -0.2756759786 -0.098116822  0.310035844
      GANIASFVMVADAMLDQGDVF            -0.6921142508 -0.233995920  0.341266313
      GCIISETGITSEQIHDIASAK            -0.2359206113 -0.045271212  0.240283188
      GGLCVDLK                                    NA           NA           NA
      ICYAFMR                          -0.4377514090 -0.058287543 -0.251418578
      NSWEGVLTGK                                  NA           NA           NA
      SLEEIVDEYSTFSESK                 -0.1884855238 -0.085688247 -0.183397329
      VLPIVSVPER                        0.2628486935  0.251265565 -0.042784340
      VTISGSGNVAQYAALK                            NA           NA           NA
      VTWENDNGEQEVAQGYR                           NA           NA           NA
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK            NA           NA           NA
      AANLGGVAVSGLEMAQNSQR             -0.1121423695  0.012272622  0.198119247
      ALVAQGVK                         -0.2531884672 -0.058050698           NA
      FIAEGSNMGSTPEAIAVFETAR           -0.1759688822 -0.053230855 -0.025955987
      GANIASFIK                        -0.1346626021  0.149556672 -0.101967636
      GCIISETGITSEQVADISSAK             0.2889078820  0.297027204           NA
      HIGQDTDVPAGDIGVGGR               -0.2395207502  0.070594907  0.192668254
      IMINCFNECIDYAK                    0.0143964776  0.078248076  0.126329445
      ITWTSER                           0.6243153621  0.677628440           NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     0.1461195038  0.399240073           NA
      SLEQIVNEYSTFSENK                  0.2343946027  0.115347817 -0.408098233
      STATGPSEAVWYGPPK                  0.0803156907  0.072060338 -0.059568315
      VDIALPCATQNEVSGEEAK               0.0022852217  0.139703709 -0.545403206
      VIELGGTVVSLSDSK                  -0.0005207847  0.134183064  0.325111208
      VQYIAGARPWTHVQK                             NA           NA           NA
      VTWENDKGEQEVAQGYR                -0.3625245313 -0.131956789  0.348974339
      AAGLTAAYAR                       -0.1351508434  0.106809508 -0.106299529
      APEAEQVLSAAATFPIAQPATDVEAR        0.2857287857 -0.360682829 -2.584903842
      AVQDNGESAFR                      -0.1302252961 -0.042380168  0.212775922
      DGKAPEAEQVLSAAATFPIAQPATDVEAR    -0.5830201862 -0.682625183 -0.531717018
      GFTLAEVK                         -0.2324508886  0.006474669  0.083150317
      IAPRPLDLLRPVVR                   -0.8586846753  0.042732320 -0.961951814
      IIVFPR                            0.8112047117  0.714976391 -0.885466541
      NQEIFDANVQR                       0.1582990486  0.168803582  0.010345525
      TIGIAVDHR                         0.4469848781  0.273836094  0.673219097
      VHFDQAGK                                    NA           NA           NA
      VHFDQAGKK                                   NA           NA           NA
      ANELLINVK                        -0.0852623550 -0.054201313  0.059594361
      ANGTTVLVGMPAGAK                             NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR          0.1983820430  0.272432269  0.424757340
      CCSDVFNQVVK                      -0.2538559864 -0.150531413  0.215911724
      DIVGAVLK                                    NA           NA           NA
      EALDFFAR                         -0.0635674360  0.053903563  0.068586586
      EKDIVGAVLK                       -0.0615652641  0.085352364  0.326554511
      GVIFYESHGK                       -0.3952059668 -0.230266752  0.821053627
      IGDYAGIK                         -0.0172961475 -0.089397052           NA
      LPLVGGHEGAGVVVGMGENVK             0.1076358647  0.248180094  0.068605265
      SANLMAGHWVAISGAAGGLGSLAVQYAK     -0.4902503314  0.305504165  0.782970287
      SIGGEVFIDFTK                     -0.2559956530 -0.121286517 -0.194051887
      SIPETQK                          -0.5481787496 -0.924604175 -0.018522232
      SISIVGSYVGNR                     -0.0141940047  0.050544896  0.016281632
      VLGIDGGEGK                        0.8316409899 -0.208213460 -0.160332082
      VLGIDGGEGKEELFR                  -0.3407958778 -0.098345213  0.287480381
      VVGLSTLPEIYEK                    -0.0260700073  0.050810083 -0.003664295
      YSGVCHTDLHAWHGDWPLPVK             0.3419087411  0.080958364 -0.033713473
      ANGTVVLVGLPAGAK                  -0.3115451482 -0.034249930  0.111166986
      CSSDVFNHVVK                                 NA           NA  0.236291866
      DIPVPKPKPNELLINVK                           NA           NA           NA
      VVGLSSLPEIYEK                     0.0771325058 -0.083091898 -0.153650898
      DIPVPEPKPNEILINVK                -0.0396578781           NA -0.250456179
      EALDFFSR                         -0.1017561194  0.165618098 -0.045093191
      GVIFYENK                                    NA           NA           NA
      IQQGTDLAEVAPILCAGVTVYK           -0.0081827733  0.071977912  0.003527630
      IVGLSELPK                        -0.1134605803  0.068536882 -0.279889231
      NMVSDIQEATK                                 NA           NA           NA
      VLGIDAGEEK                                  NA           NA           NA
                                       logRatio_3_4 logRatio_3_5 logRatio_3_6
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                  -0.142439760   0.61616938 -0.322140001
      AEWALR                           -0.152991860  -0.03981807 -0.026239789
      DEGLHTDFACLLFAHLK                -0.120617530           NA           NA
      DIHDWNNR                                   NA           NA           NA
      ELETLREENR                                 NA           NA           NA
      ESEFLFNAIHTIPEIGEK               -0.107301956  -0.22113043  0.160007591
      GMMPGLTFSNELICR                   0.085091781           NA -0.025250920
      IVTEAVEIEQR                      -0.342898412  -0.33608002  0.225697655
      LLVAFGNK                         -0.407012842  -0.98529653  0.313518981
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                       -0.246695322  -0.02048541  0.419447324
      TNFFEK                           -0.194549576  -0.25241139 -0.172513042
      TVLFPIK                          -0.352774062  -0.49454306  0.200870290
      VENPFDFMENISLAGK                 -0.119590999   0.01970767 -0.246935371
      WIQDADALFGER                      0.032321287  -0.33250274 -0.241573991
      YFLDALPVALLGMNADLMNQYVEFVADR     -0.150861222  -0.18460990  0.019210276
      AANLGGVAVSGLEMAQNSQK             -0.022548555  -0.15306288  0.038215791
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                      -0.237230938  -0.33419126  0.473949377
      FLGFEQIFK                         0.045080743   0.20915400  0.156012000
      GANIASFVMVADAMLDQGDVF             0.123270962   0.65091763 -0.063771312
      GCIISETGITSEQIHDIASAK            -0.128489834  -0.01174687  0.152239263
      GGLCVDLK                         -0.250316929  -0.10969772  0.172441571
      ICYAFMR                           0.056973769   0.36158670  0.209507050
      NSWEGVLTGK                       -0.260453691  -0.16512543  0.129695374
      SLEEIVDEYSTFSESK                 -0.005922440  -0.13695467 -0.440924685
      VLPIVSVPER                       -0.200145586  -0.41924544 -0.031997566
      VTISGSGNVAQYAALK                           NA           NA           NA
      VTWENDNGEQEVAQGYR                          NA           NA           NA
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR              0.018897356   0.18283349 -0.153998296
      ALVAQGVK                                   NA  -0.39151144 -0.027520730
      FIAEGSNMGSTPEAIAVFETAR           -0.219659522  -0.18684484 -0.485433581
      GANIASFIK                                  NA           NA           NA
      GCIISETGITSEQVADISSAK             0.048604186   0.15602879  0.436771245
      HIGQDTDVPAGDIGVGGR               -0.182046720  -0.08779136  0.371912852
      IMINCFNECIDYAK                    0.117409519   0.13969561  0.176296078
      ITWTSER                                    NA  -0.35173187           NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR    -0.002104056   0.15906890 -1.314892175
      SLEQIVNEYSTFSENK                 -0.088950144  -0.06333726  0.004788589
      STATGPSEAVWYGPPK                 -0.005646720  -0.15563822           NA
      VDIALPCATQNEVSGEEAK              -0.007101587  -0.21818907  0.024789660
      VIELGGTVVSLSDSK                  -0.206740609  -0.26298061 -0.089775556
      VQYIAGARPWTHVQK                            NA   0.84739503           NA
      VTWENDKGEQEVAQGYR                -0.037602660   0.19465260  0.296495484
      AAGLTAAYAR                       -0.336208150  -0.23454246  0.080684482
      APEAEQVLSAAATFPIAQPATDVEAR       -0.278730248  -0.92384953 -0.070446209
      AVQDNGESAFR                      -0.115147280  -0.72149400  0.116150235
      DGKAPEAEQVLSAAATFPIAQPATDVEAR    -0.003534379   0.20852487 -0.659973222
      GFTLAEVK                         -0.013483010   0.05568670  0.182000368
      IAPRPLDLLRPVVR                   -0.486812450   0.19709139 -0.159762640
      IIVFPR                           -0.072857191  -0.38084081  0.068002740
      NQEIFDANVQR                      -0.104963902  -0.16121876  0.054885817
      TIGIAVDHR                        -0.316175071  -0.39751927           NA
      VHFDQAGK                                   NA  -0.22675374  0.936901091
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                        -0.366573686  -0.25483720  0.317266877
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR         -0.275197129  -0.47915640  0.472640216
      CCSDVFNQVVK                       0.001001409   0.06273752  0.138231890
      DIVGAVLK                                   NA           NA  0.149565914
      EALDFFAR                         -0.090597634  -0.12119455  0.056415760
      EKDIVGAVLK                       -0.168492098  -0.01072552  0.368190336
      GVIFYESHGK                       -0.214111923  -0.07526648  0.448330193
      IGDYAGIK                         -0.362471646           NA -0.040954654
      LPLVGGHEGAGVVVGMGENVK             0.046886153  -0.13936817  0.343589651
      SANLMAGHWVAISGAAGGLGSLAVQYAK     -0.380321329  -0.39618426 -0.464029733
      SIGGEVFIDFTK                     -0.099333002  -0.24452581 -0.194505291
      SIPETQK                           0.426902334   0.26475902 -0.804582609
      SISIVGSYVGNR                     -0.199213653  -0.24712181  0.171028780
      VLGIDGGEGK                       -0.522551645  -1.44795414 -0.025635028
      VLGIDGGEGKEELFR                   0.016931703   0.12596301  0.269279607
      VVGLSTLPEIYEK                    -0.103936949  -0.10335674  0.007150177
      YSGVCHTDLHAWHGDWPLPVK             0.424950098  -0.03191982  0.071340523
      ANGTVVLVGLPAGAK                  -0.122148889  -0.04863585  0.214877974
      CSSDVFNHVVK                                NA           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                    -0.247729760  -0.31205006 -0.112419817
      DIPVPEPKPNEILINVK                -0.120482589  -0.80648170 -0.141539191
      EALDFFSR                          0.181057126           NA           NA
      GVIFYENK                         -0.026019431           NA           NA
      IQQGTDLAEVAPILCAGVTVYK           -0.061301180   0.13146558 -0.037649381
      IVGLSELPK                         0.020511007  -0.03533411  0.241181439
      NMVSDIQEATK                      -0.065464841  -0.16798770 -0.188908812
      VLGIDAGEEK                                 NA           NA           NA
                                       logRatio_3_7 logRatio_3_8 logRatio_3_9
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                            NA -0.213460040           NA
      AEWALR                           -0.084220443 -0.013211444  0.058516572
      DEGLHTDFACLLFAHLK                          NA           NA  0.938116572
      DIHDWNNR                                   NA           NA           NA
      ELETLREENR                                 NA           NA           NA
      ESEFLFNAIHTIPEIGEK               -0.113401397 -0.006913222  0.006728633
      GMMPGLTFSNELICR                   0.263276882  0.108804639           NA
      IVTEAVEIEQR                      -0.425651996 -0.298039580  0.198125202
      LLVAFGNK                                   NA -0.632405387 -0.822843336
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                       -0.160620484  0.141651714  0.171184492
      TNFFEK                           -0.084038772 -0.138531059 -0.110994701
      TVLFPIK                          -0.329109724 -0.399157843 -0.393837940
      VENPFDFMENISLAGK                  0.028623441 -0.033586092 -0.210233143
      WIQDADALFGER                      0.045296709 -0.253003814 -0.295284378
      YFLDALPVALLGMNADLMNQYVEFVADR     -0.217105091 -0.032463778  0.217720472
      AANLGGVAVSGLEMAQNSQK             -0.289322802  0.082534414           NA
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                      -0.471440690 -0.468874077  0.558900122
      FLGFEQIFK                         0.069365611  0.246924767  0.655077433
      GANIASFVMVADAMLDQGDVF             0.073645777  0.531764107  1.107026340
      GCIISETGITSEQIHDIASAK            -0.176305886  0.014343514  0.299897914
      GGLCVDLK                                   NA           NA -0.016219476
      ICYAFMR                           0.138176891  0.517640756  0.324509722
      NSWEGVLTGK                       -0.349967919 -0.099488192 -0.065634720
      SLEEIVDEYSTFSESK                 -0.091694673  0.011102604 -0.086606478
      VLPIVSVPER                       -0.108567048 -0.120150176 -0.414200082
      VTISGSGNVAQYAALK                           NA           NA           NA
      VTWENDNGEQEVAQGYR                          NA           NA           NA
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR             -0.357663746 -0.233248755 -0.047402130
      ALVAQGVK                         -0.455274857 -0.260137089           NA
      FIAEGSNMGSTPEAIAVFETAR           -0.121522766  0.001215261  0.028490129
      GANIASFIK                                  NA           NA           NA
      GCIISETGITSEQVADISSAK             0.145416474  0.153535796           NA
      HIGQDTDVPAGDIGVGGR               -0.225493865  0.084621792  0.206695140
      IMINCFNECIDYAK                    0.110835075  0.174686673  0.222768043
      ITWTSER                          -0.188307562 -0.134994484           NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR    -0.082482984  0.170637585           NA
      SLEQIVNEYSTFSENK                  0.014923106 -0.104123679 -0.627569729
      STATGPSEAVWYGPPK                 -0.017279918 -0.025535271 -0.157163923
      VDIALPCATQNEVSGEEAK              -0.071933484  0.065485003 -0.619621912
      VIELGGTVVSLSDSK                  -0.216327577 -0.081623729  0.109304415
      VQYIAGARPWTHVQK                   0.980264121           NA           NA
      VTWENDKGEQEVAQGYR                 0.036223893  0.266791636  0.747722764
      AAGLTAAYAR                       -0.278337712 -0.036377360 -0.249486397
      APEAEQVLSAAATFPIAQPATDVEAR       -0.273289293 -0.919700908 -3.143921921
      AVQDNGESAFR                      -0.078425903  0.009419225  0.264575315
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     0.277921830  0.178316833  0.329224999
      GFTLAEVK                         -0.082509728  0.156415830  0.233091478
      IAPRPLDLLRPVVR                   -0.234625497  0.666791498 -0.337892636
      IIVFPR                           -0.098589238 -0.194817559 -1.795260491
      NQEIFDANVQR                       0.001916750  0.012421283 -0.146036774
      TIGIAVDHR                        -0.069012065 -0.242160849  0.157222154
      VHFDQAGK                         -0.413850332 -0.115905092  0.284160983
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                        -0.240408660 -0.209347618 -0.095551944
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR         -0.178674376 -0.104624150  0.047700922
      CCSDVFNQVVK                      -0.067274381  0.036050193  0.402493330
      DIVGAVLK                                   NA -0.062581855  0.085070481
      EALDFFAR                         -0.076293192  0.041177806  0.055860830
      EKDIVGAVLK                       -0.109832352  0.037085276  0.278287422
      GVIFYESHGK                       -0.174199917 -0.009260702  1.042059677
      IGDYAGIK                         -0.330650323 -0.402751228           NA
      LPLVGGHEGAGVVVGMGENVK             0.107183462  0.247727692  0.068152862
      SANLMAGHWVAISGAAGGLGSLAVQYAK     -0.860445606 -0.064691110  0.412775012
      SIGGEVFIDFTK                     -0.138630758 -0.003921622 -0.076686992
      SIPETQK                           0.443671146  0.067245720  0.973327663
      SISIVGSYVGNR                     -0.124049419 -0.059310518 -0.093573782
      VLGIDGGEGK                       -0.503030828 -1.542885279 -1.495003900
      VLGIDGGEGKEELFR                   0.009620077  0.252070742  0.637896335
      VVGLSTLPEIYEK                    -0.030160560  0.046719530 -0.007754847
      YSGVCHTDLHAWHGDWPLPVK             0.322395582  0.061445205 -0.053226633
      ANGTVVLVGLPAGAK                  -0.205991093  0.071304125  0.216721041
      CSSDVFNHVVK                                NA           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                    -0.118108178 -0.278332581 -0.348891581
      DIPVPEPKPNEILINVK                -0.254123786           NA -0.464922087
      EALDFFSR                          0.021637509  0.289011726  0.078300437
      GVIFYENK                                   NA  0.064811646  0.425680552
      IQQGTDLAEVAPILCAGVTVYK            0.095266074  0.175426759  0.106976476
      IVGLSELPK                        -0.003133325  0.178864137 -0.169561976
      NMVSDIQEATK                      -0.124262470 -0.134199571  0.187771055
      VLGIDAGEEK                                 NA           NA           NA
                                        logRatio_4_5 logRatio_4_6 logRatio_4_7
      AAADALSDLEIK                                NA           NA           NA
      AAADALSDLEIKDSK                   0.7586091393  -0.17970024           NA
      AEWALR                            0.1131737882   0.12675207  0.068771417
      DEGLHTDFACLLFAHLK                           NA           NA           NA
      DIHDWNNR                                    NA           NA           NA
      ELETLREENR                                  NA           NA           NA
      ESEFLFNAIHTIPEIGEK               -0.1138284719   0.26730955 -0.006099441
      GMMPGLTFSNELICR                             NA  -0.11034270  0.178185101
      IVTEAVEIEQR                       0.0068183919   0.56859607 -0.082753584
      LLVAFGNK                         -0.5782836843   0.72053182           NA
      LLVAFGNKK                                   NA           NA           NA
      NKPDPAIVEK                        0.2262099120   0.66614265  0.086074838
      TNFFEK                           -0.0578618141   0.02203653  0.110510804
      TVLFPIK                          -0.1417689950   0.55364435  0.023664338
      VENPFDFMENISLAGK                  0.1392986706  -0.12734437  0.148214440
      WIQDADALFGER                     -0.3648240286  -0.27389528  0.012975422
      YFLDALPVALLGMNADLMNQYVEFVADR     -0.0337486749   0.17007150 -0.066243868
      AANLGGVAVSGLEMAQNSQK             -0.1305143247   0.06076435 -0.266774247
      DAVWFGPPK                                   NA           NA           NA
      EIGYLFGAYR                                  NA           NA           NA
      FHPSVNLSILK                      -0.0969603182   0.71118031 -0.234209753
      FLGFEQIFK                         0.1640732599   0.11093126  0.024284868
      GANIASFVMVADAMLDQGDVF             0.5276466700  -0.18704227 -0.049625185
      GCIISETGITSEQIHDIASAK             0.1167429660   0.28072910 -0.047816052
      GGLCVDLK                          0.1406192098   0.42275850           NA
      ICYAFMR                           0.3046129354   0.15253328  0.081203122
      NSWEGVLTGK                        0.0953282562   0.39014906 -0.089514228
      SLEEIVDEYSTFSESK                 -0.1310322259  -0.43500225 -0.085772233
      VLPIVSVPER                       -0.2190998588   0.16814802  0.091578538
      VTISGSGNVAQYAALK                            NA           NA           NA
      VTWENDNGEQEVAQGYR                -0.0952666731           NA  0.170226342
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK            NA           NA           NA
      AANLGGVAVSGLEMAQNSQR              0.1639361325  -0.17289565 -0.376561103
      ALVAQGVK                                    NA           NA           NA
      FIAEGSNMGSTPEAIAVFETAR            0.0328146790  -0.26577406  0.098136756
      GANIASFIK                         0.0275088712   0.26527384  0.011721746
      GCIISETGITSEQVADISSAK             0.1074246072   0.38816706  0.096812287
      HIGQDTDVPAGDIGVGGR                0.0942553630   0.55395957 -0.043447145
      IMINCFNECIDYAK                    0.0222860913   0.05888656 -0.006574445
      ITWTSER                                     NA           NA           NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     0.1611729584  -1.31278812 -0.080378928
      SLEQIVNEYSTFSENK                  0.0256128839   0.09373873  0.103873250
      STATGPSEAVWYGPPK                 -0.1499915021           NA -0.011633198
      VDIALPCATQNEVSGEEAK              -0.2110874851   0.03189125 -0.064831897
      VIELGGTVVSLSDSK                  -0.0562399965   0.11696505 -0.009586969
      VQYIAGARPWTHVQK                             NA           NA           NA
      VTWENDKGEQEVAQGYR                 0.2322552597   0.33409814  0.073826553
      AAGLTAAYAR                        0.1016656948   0.41689263  0.057870439
      APEAEQVLSAAATFPIAQPATDVEAR       -0.6451192797   0.20828404  0.005440955
      AVQDNGESAFR                      -0.6063467160   0.23129752  0.036721377
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     0.2120592509  -0.65643884  0.281456210
      GFTLAEVK                          0.0691697058   0.19548338 -0.069026717
      IAPRPLDLLRPVVR                    0.6839038435   0.32704981  0.252186953
      IIVFPR                           -0.3079836224   0.14085993 -0.025732047
      NQEIFDANVQR                      -0.0562548567   0.15984972  0.106880652
      TIGIAVDHR                        -0.0813442022           NA  0.247163005
      VHFDQAGK                                    NA           NA           NA
      VHFDQAGKK                                   NA           NA           NA
      ANELLINVK                         0.1117364892   0.68384056  0.126165026
      ANGTTVLVGMPAGAK                             NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR         -0.2039592741   0.74783735  0.096522753
      CCSDVFNQVVK                       0.0617361062   0.13723048 -0.068275790
      DIVGAVLK                                    NA           NA           NA
      EALDFFAR                         -0.0305969206   0.14701339  0.014304442
      EKDIVGAVLK                        0.1577665827   0.53668243  0.058659745
      GVIFYESHGK                        0.1388454469   0.66244212  0.039912006
      IGDYAGIK                                    NA   0.32151699  0.031821323
      LPLVGGHEGAGVVVGMGENVK            -0.1862543251   0.29670350  0.060297309
      SANLMAGHWVAISGAAGGLGSLAVQYAK     -0.0158629303  -0.08370840 -0.480124277
      SIGGEVFIDFTK                     -0.1451928065  -0.09517229 -0.039297756
      SIPETQK                          -0.1621433112  -1.23148494  0.016768811
      SISIVGSYVGNR                     -0.0479081539   0.37024243  0.075164234
      VLGIDGGEGK                       -0.9254024986   0.49691662  0.019520817
      VLGIDGGEGKEELFR                   0.1090313032   0.25234790 -0.007311627
      VVGLSTLPEIYEK                     0.0005802112   0.11108713  0.073776389
      YSGVCHTDLHAWHGDWPLPVK            -0.4568699154  -0.35360958 -0.102554516
      ANGTVVLVGLPAGAK                   0.0735130362   0.33702686 -0.083842204
      CSSDVFNHVVK                                 NA           NA           NA
      DIPVPKPKPNELLINVK                           NA           NA           NA
      VVGLSSLPEIYEK                    -0.0643202998   0.13530994  0.129621583
      DIPVPEPKPNEILINVK                -0.6859991144  -0.02105660 -0.133641197
      EALDFFSR                                    NA           NA -0.159419617
      GVIFYENK                                    NA           NA           NA
      IQQGTDLAEVAPILCAGVTVYK            0.1927667565   0.02365180  0.156567253
      IVGLSELPK                        -0.0558451144   0.22067043 -0.023644332
      NMVSDIQEATK                      -0.1025228571  -0.12344397 -0.058797629
      VLGIDAGEEK                                  NA           NA           NA
                                       logRatio_4_8 logRatio_4_9 logRatio_5_6
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                   -0.07102028           NA  -0.93830938
      AEWALR                             0.13978042   0.21150843   0.01357828
      DEGLHTDFACLLFAHLK                          NA   1.05873410           NA
      DIHDWNNR                                   NA           NA           NA
      ELETLREENR                                 NA           NA   0.52754785
      ESEFLFNAIHTIPEIGEK                 0.10038873   0.11403059   0.38113802
      GMMPGLTFSNELICR                    0.02371286           NA           NA
      IVTEAVEIEQR                        0.04485883   0.54102361   0.56177767
      LLVAFGNK                          -0.22539255  -0.41583049   1.29881551
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                         0.38834704   0.41787981   0.43993273
      TNFFEK                             0.05601852   0.08355488   0.07989835
      TVLFPIK                           -0.04638378  -0.04106388   0.69541335
      VENPFDFMENISLAGK                   0.08600491  -0.09064214  -0.26664304
      WIQDADALFGER                      -0.28532510  -0.32760566   0.09092875
      YFLDALPVALLGMNADLMNQYVEFVADR       0.11839744   0.36858169   0.20382017
      AANLGGVAVSGLEMAQNSQK               0.10508297           NA   0.19127867
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                       -0.23164314   0.79613106   0.80814063
      FLGFEQIFK                          0.20184402   0.60999669  -0.05314200
      GANIASFVMVADAMLDQGDVF              0.40849314   0.98375538  -0.71468894
      GCIISETGITSEQIHDIASAK              0.14283335   0.42838775   0.16398613
      GGLCVDLK                                   NA   0.23409745   0.28213929
      ICYAFMR                            0.46066699   0.26753595  -0.15207965
      NSWEGVLTGK                         0.16096550   0.19481897   0.29482081
      SLEEIVDEYSTFSESK                   0.01702504  -0.08068404  -0.30397002
      VLPIVSVPER                         0.07999541  -0.21405450   0.38724788
      VTISGSGNVAQYAALK                           NA           NA           NA
      VTWENDNGEQEVAQGYR                 -0.11366378   0.10916431           NA
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR              -0.25214611  -0.06629949  -0.33683178
      ALVAQGVK                                   NA           NA   0.36399071
      FIAEGSNMGSTPEAIAVFETAR             0.22087478   0.24814965  -0.29858874
      GANIASFIK                          0.29594102   0.04441671   0.23776497
      GCIISETGITSEQVADISSAK              0.10493161           NA   0.28074245
      HIGQDTDVPAGDIGVGGR                 0.26666851   0.38874186   0.45970421
      IMINCFNECIDYAK                     0.05727715   0.10535852   0.03660047
      ITWTSER                                    NA           NA           NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR      0.17274164           NA  -1.47396108
      SLEQIVNEYSTFSENK                  -0.01517353  -0.53861958   0.06812585
      STATGPSEAVWYGPPK                  -0.01988855  -0.15151720           NA
      VDIALPCATQNEVSGEEAK                0.07258659  -0.61252032   0.24297873
      VIELGGTVVSLSDSK                    0.12511688   0.31604502   0.17320505
      VQYIAGARPWTHVQK                            NA           NA           NA
      VTWENDKGEQEVAQGYR                  0.30439430   0.78532542   0.10184288
      AAGLTAAYAR                         0.29983079   0.08672175   0.31522694
      APEAEQVLSAAATFPIAQPATDVEAR        -0.64097066  -2.86519167   0.85340332
      AVQDNGESAFR                        0.12456650   0.37972259   0.83764423
      DGKAPEAEQVLSAAATFPIAQPATDVEAR      0.18185121   0.33275938  -0.86849809
      GFTLAEVK                           0.16989884   0.24657449   0.12631367
      IAPRPLDLLRPVVR                     1.15360395   0.14891981  -0.35685403
      IIVFPR                            -0.12196037  -1.72240330   0.44884355
      NQEIFDANVQR                        0.11738518  -0.04107287   0.21610458
      TIGIAVDHR                          0.07401422   0.47339722           NA
      VHFDQAGK                                   NA           NA   1.16365484
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                          0.15722607   0.27102174   0.57210407
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR           0.17057298   0.32289805   0.95179662
      CCSDVFNQVVK                        0.03504878   0.40149192   0.07549437
      DIVGAVLK                                   NA           NA           NA
      EALDFFAR                           0.13177544   0.14645846   0.17761031
      EKDIVGAVLK                         0.20557737   0.44677952   0.37891585
      GVIFYESHGK                         0.20485122   1.25617160   0.52359667
      IGDYAGIK                          -0.04027958           NA           NA
      LPLVGGHEGAGVVVGMGENVK              0.20084154   0.02126671   0.48295782
      SANLMAGHWVAISGAAGGLGSLAVQYAK       0.31563022   0.79309634  -0.06784547
      SIGGEVFIDFTK                       0.09541138   0.02264601   0.05002052
      SIPETQK                           -0.35965661   0.54642533  -1.06934163
      SISIVGSYVGNR                       0.13990313   0.10563987   0.41815059
      VLGIDGGEGK                        -1.02033363  -0.97245225   1.42231912
      VLGIDGGEGKEELFR                    0.23513904   0.62096463   0.14331660
      VVGLSTLPEIYEK                      0.15065648   0.09618210   0.11050691
      YSGVCHTDLHAWHGDWPLPVK             -0.36350489  -0.47817673   0.10326034
      ANGTVVLVGLPAGAK                    0.19345301   0.33886993   0.26351383
      CSSDVFNHVVK                                NA           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                     -0.03060282  -0.10116182   0.19963024
      DIPVPEPKPNEILINVK                          NA  -0.34443950   0.66494251
      EALDFFSR                           0.10795460  -0.10275669           NA
      GVIFYENK                           0.09083108   0.45169998           NA
      IQQGTDLAEVAPILCAGVTVYK             0.23672794   0.16827766  -0.16911496
      IVGLSELPK                          0.15835313  -0.19007298   0.27651555
      NMVSDIQEATK                       -0.06873473   0.25323590  -0.02092111
      VLGIDAGEEK                                 NA           NA           NA
                                       logRatio_5_7 logRatio_5_8 logRatio_5_9
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                            NA -0.829629419           NA
      AEWALR                           -0.044402371  0.026606629  0.098334644
      DEGLHTDFACLLFAHLK                          NA           NA           NA
      DIHDWNNR                                   NA           NA           NA
      ELETLREENR                                 NA  0.215235341  1.020989595
      ESEFLFNAIHTIPEIGEK                0.107729031  0.214217206  0.227859061
      GMMPGLTFSNELICR                            NA           NA           NA
      IVTEAVEIEQR                      -0.089571976  0.038040439  0.534205222
      LLVAFGNK                                   NA  0.352891139  0.162453190
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                       -0.140135074  0.162137124  0.191669902
      TNFFEK                            0.168372618  0.113880331  0.141416689
      TVLFPIK                           0.165433333  0.095385214  0.100705117
      VENPFDFMENISLAGK                  0.008915769 -0.053293763 -0.229940814
      WIQDADALFGER                      0.377799451  0.079498927  0.037218364
      YFLDALPVALLGMNADLMNQYVEFVADR     -0.032495194  0.152146119  0.402330370
      AANLGGVAVSGLEMAQNSQK             -0.136259922  0.235597294           NA
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                      -0.137249435 -0.134682821  0.893091378
      FLGFEQIFK                        -0.139788392  0.037770764  0.445923431
      GANIASFVMVADAMLDQGDVF            -0.577271855 -0.119153525  0.456108708
      GCIISETGITSEQIHDIASAK            -0.164559018  0.026090382  0.311644782
      GGLCVDLK                                   NA           NA  0.093478243
      ICYAFMR                          -0.223409814  0.156054052 -0.037076983
      NSWEGVLTGK                       -0.184842484  0.065637243  0.099490715
      SLEEIVDEYSTFSESK                  0.045259993  0.148057270  0.050348188
      VLPIVSVPER                        0.310678397  0.299095268  0.005045363
      VTISGSGNVAQYAALK                           NA           NA           NA
      VTWENDNGEQEVAQGYR                 0.265493015 -0.018397103  0.204430980
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR             -0.540497235 -0.416082244 -0.230235619
      ALVAQGVK                         -0.063763414  0.131374355           NA
      FIAEGSNMGSTPEAIAVFETAR            0.065322077  0.188060104  0.215334972
      GANIASFIK                        -0.015787125  0.268432149  0.016907841
      GCIISETGITSEQVADISSAK            -0.010612320 -0.002492998           NA
      HIGQDTDVPAGDIGVGGR               -0.137702508  0.172413149  0.294486496
      IMINCFNECIDYAK                   -0.028860536  0.034991062  0.083072432
      ITWTSER                           0.163424309  0.216737387           NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR    -0.241551886  0.011568683           NA
      SLEQIVNEYSTFSENK                  0.078260367 -0.040786419 -0.564232469
      STATGPSEAVWYGPPK                  0.138358304  0.130102951 -0.001525701
      VDIALPCATQNEVSGEEAK               0.146255589  0.283674076 -0.401432839
      VIELGGTVVSLSDSK                   0.046653028  0.181356876  0.372285020
      VQYIAGARPWTHVQK                   0.132869095           NA           NA
      VTWENDKGEQEVAQGYR                -0.158428707  0.072139036  0.553070164
      AAGLTAAYAR                       -0.043795256  0.198165095 -0.014943942
      APEAEQVLSAAATFPIAQPATDVEAR        0.650560235  0.004148620 -2.220072393
      AVQDNGESAFR                       0.643068093  0.730913221  0.986069311
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     0.069396959 -0.030208038  0.120700127
      GFTLAEVK                         -0.138196423  0.100729134  0.177404783
      IAPRPLDLLRPVVR                   -0.431716891  0.469700104 -0.534984029
      IIVFPR                            0.282251575  0.186023254 -1.414419678
      NQEIFDANVQR                       0.163135509  0.173640042  0.015181985
      TIGIAVDHR                         0.328507208  0.155358424  0.554741427
      VHFDQAGK                         -0.187096588  0.110848652  0.510914727
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                         0.014428537  0.045489579  0.159285253
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR          0.300482027  0.374532253  0.526857325
      CCSDVFNQVVK                      -0.130011896 -0.026687322  0.339755814
      DIVGAVLK                                   NA           NA           NA
      EALDFFAR                          0.044901362  0.162372361  0.177055385
      EKDIVGAVLK                       -0.099106837  0.047810791  0.289012937
      GVIFYESHGK                       -0.098933441  0.066005774  1.117326153
      IGDYAGIK                                   NA           NA           NA
      LPLVGGHEGAGVVVGMGENVK             0.246551635  0.387095864  0.207521035
      SANLMAGHWVAISGAAGGLGSLAVQYAK     -0.464261347  0.331493149  0.808959271
      SIGGEVFIDFTK                      0.105895051  0.240604187  0.167838817
      SIPETQK                           0.178912122 -0.197513303  0.708568640
      SISIVGSYVGNR                      0.123072388  0.187811289  0.153548024
      VLGIDGGEGK                        0.944923315 -0.094931135 -0.047049756
      VLGIDGGEGKEELFR                  -0.116342930  0.126107735  0.511933329
      VVGLSTLPEIYEK                     0.073196178  0.150076268  0.095601890
      YSGVCHTDLHAWHGDWPLPVK             0.354315399  0.093365022 -0.021306815
      ANGTVVLVGLPAGAK                  -0.157355240  0.119939978  0.265356894
      CSSDVFNHVVK                                NA           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                     0.193941883  0.033717479 -0.036841521
      DIPVPEPKPNEILINVK                 0.552357917           NA  0.341559616
      EALDFFSR                                   NA           NA           NA
      GVIFYENK                                   NA           NA           NA
      IQQGTDLAEVAPILCAGVTVYK           -0.036199503  0.043961182 -0.024489100
      IVGLSELPK                         0.032200783  0.214198245 -0.134227868
      NMVSDIQEATK                       0.043725228  0.033788127  0.355758753
      VLGIDAGEEK                                 NA           NA           NA
                                       logRatio_6_7 logRatio_6_8  logRatio_6_9
      AAADALSDLEIK                               NA           NA            NA
      AAADALSDLEIKDSK                            NA  0.108679962            NA
      AEWALR                           -0.057980654  0.013028346  0.0847563614
      DEGLHTDFACLLFAHLK                          NA           NA            NA
      DIHDWNNR                                   NA           NA            NA
      ELETLREENR                                 NA -0.312312508  0.4934417459
      ESEFLFNAIHTIPEIGEK               -0.273408987 -0.166920813 -0.1532789576
      GMMPGLTFSNELICR                   0.288527802  0.134055558            NA
      IVTEAVEIEQR                      -0.651349651 -0.523737235 -0.0275724526
      LLVAFGNK                                   NA -0.945924368 -1.1363623169
      LLVAFGNKK                                  NA           NA            NA
      NKPDPAIVEK                       -0.580067808 -0.277795609 -0.2482628315
      TNFFEK                            0.088474270  0.033981983  0.0615183415
      TVLFPIK                          -0.529980015 -0.600028133 -0.5947082302
      VENPFDFMENISLAGK                  0.275558812  0.213349279  0.0367022279
      WIQDADALFGER                      0.286870700 -0.011429824 -0.0537103870
      YFLDALPVALLGMNADLMNQYVEFVADR     -0.236315366 -0.051674054  0.1985101967
      AANLGGVAVSGLEMAQNSQK             -0.327538592  0.044318623            NA
      DAVWFGPPK                                  NA           NA            NA
      EIGYLFGAYR                                 NA           NA            NA
      FHPSVNLSILK                      -0.945390067 -0.942823454  0.0849507448
      FLGFEQIFK                        -0.086646389  0.090912767  0.4990654338
      GANIASFVMVADAMLDQGDVF             0.137417089  0.595535419  1.1707976525
      GCIISETGITSEQIHDIASAK            -0.328545148 -0.137895749  0.1476586510
      GGLCVDLK                                   NA           NA -0.1886610468
      ICYAFMR                          -0.071330160  0.308133706  0.1150026714
      NSWEGVLTGK                       -0.479663293 -0.229183566 -0.1953300939
      SLEEIVDEYSTFSESK                  0.349230012  0.452027289  0.3543182069
      VLPIVSVPER                       -0.076569481 -0.088152610 -0.3822025153
      VTISGSGNVAQYAALK                           NA           NA            NA
      VTWENDNGEQEVAQGYR                          NA           NA            NA
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA            NA
      AANLGGVAVSGLEMAQNSQR             -0.203665451 -0.079250460  0.1065961658
      ALVAQGVK                         -0.427754128 -0.232616359            NA
      FIAEGSNMGSTPEAIAVFETAR            0.363910816  0.486648843  0.5139237107
      GANIASFIK                        -0.253552095  0.030667179 -0.2208571288
      GCIISETGITSEQVADISSAK            -0.291354772 -0.283235450            NA
      HIGQDTDVPAGDIGVGGR               -0.597406717 -0.287291060 -0.1652177126
      IMINCFNECIDYAK                   -0.065461004 -0.001609406  0.0464719641
      ITWTSER                                    NA           NA            NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     1.232409192  1.485529761            NA
      SLEQIVNEYSTFSENK                  0.010134517 -0.108912268 -0.6323583184
      STATGPSEAVWYGPPK                           NA           NA            NA
      VDIALPCATQNEVSGEEAK              -0.096723144  0.040695343 -0.6444115719
      VIELGGTVVSLSDSK                  -0.126552022  0.008151827  0.1990799708
      VQYIAGARPWTHVQK                            NA           NA            NA
      VTWENDKGEQEVAQGYR                -0.260271591 -0.029703848  0.4512272795
      AAGLTAAYAR                       -0.359022194 -0.117061842 -0.3301708789
      APEAEQVLSAAATFPIAQPATDVEAR       -0.202843084 -0.849254699 -3.0734757119
      AVQDNGESAFR                      -0.194576138 -0.106731010  0.1484250795
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     0.937895053  0.838290056  0.9891982207
      GFTLAEVK                         -0.264510096 -0.025584538  0.0510911104
      IAPRPLDLLRPVVR                   -0.074862856  0.826554138 -0.1781299952
      IIVFPR                           -0.166591978 -0.262820299 -1.8632632309
      NQEIFDANVQR                      -0.052969067 -0.042464534 -0.2009225905
      TIGIAVDHR                                  NA           NA            NA
      VHFDQAGK                         -1.350751424 -1.052806184 -0.6527401085
      VHFDQAGKK                                  NA           NA            NA
      ANELLINVK                        -0.557675537 -0.526614495 -0.4128188213
      ANGTTVLVGMPAGAK                  -0.157713214           NA -0.0869716916
      ATDGGAHGVINVSVSEAAIEASTR         -0.651314592 -0.577264367 -0.4249392948
      CCSDVFNQVVK                      -0.205506271 -0.102181697  0.2642614394
      DIVGAVLK                                   NA -0.212147769 -0.0644954327
      EALDFFAR                         -0.132708952 -0.015237953 -0.0005549297
      EKDIVGAVLK                       -0.478022688 -0.331105060 -0.0899029134
      GVIFYESHGK                       -0.622530109 -0.457590894  0.5937294846
      IGDYAGIK                         -0.289695669 -0.361796574            NA
      LPLVGGHEGAGVVVGMGENVK            -0.236406188 -0.095861959 -0.2754367881
      SANLMAGHWVAISGAAGGLGSLAVQYAK     -0.396415873  0.399338623  0.8768047450
      SIGGEVFIDFTK                      0.055874533  0.190583669  0.1178182988
      SIPETQK                           1.248253754  0.871828329  1.7779102719
      SISIVGSYVGNR                     -0.295078199 -0.230339298 -0.2646025621
      VLGIDGGEGK                       -0.477395800 -1.517250250 -1.4693688713
      VLGIDGGEGKEELFR                  -0.259659530 -0.017208865  0.3686167283
      VVGLSTLPEIYEK                    -0.037310736  0.039569354 -0.0149050239
      YSGVCHTDLHAWHGDWPLPVK             0.251055059 -0.009895318 -0.1245671556
      ANGTVVLVGLPAGAK                  -0.420869067 -0.143573849  0.0018430670
      CSSDVFNHVVK                                NA           NA  0.2965473570
      DIPVPKPKPNELLINVK                          NA           NA            NA
      VVGLSSLPEIYEK                    -0.005688361 -0.165912764 -0.2364717643
      DIPVPEPKPNEILINVK                -0.112584596           NA -0.3233828963
      EALDFFSR                                   NA           NA            NA
      GVIFYENK                                   NA           NA            NA
      IQQGTDLAEVAPILCAGVTVYK            0.132915454  0.213076139  0.1446258571
      IVGLSELPK                        -0.244314764 -0.062317302 -0.4107434146
      NMVSDIQEATK                       0.064646342  0.054709240  0.3766798668
      VLGIDAGEEK                                 NA           NA            NA
                                       logRatio_7_8 logRatio_7_9 logRatio_8_9
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                            NA           NA           NA
      AEWALR                            0.071008999  0.142737015  0.071728016
      DEGLHTDFACLLFAHLK                          NA           NA           NA
      DIHDWNNR                                   NA           NA           NA
      ELETLREENR                                 NA           NA  0.805754254
      ESEFLFNAIHTIPEIGEK                0.106488175  0.120130030  0.013641855
      GMMPGLTFSNELICR                  -0.154472243           NA           NA
      IVTEAVEIEQR                       0.127612415  0.623777198  0.496164783
      LLVAFGNK                                   NA           NA -0.190437948
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                        0.302272198  0.331804976  0.029532778
      TNFFEK                           -0.054492287 -0.026955929  0.027536358
      TVLFPIK                          -0.070048119 -0.064728216  0.005319903
      VENPFDFMENISLAGK                 -0.062209533 -0.238856584 -0.176647051
      WIQDADALFGER                     -0.298300524 -0.340581087 -0.042280563
      YFLDALPVALLGMNADLMNQYVEFVADR      0.184641312  0.434825563  0.250184251
      AANLGGVAVSGLEMAQNSQK              0.371857216           NA           NA
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                       0.002566614  1.030340812  1.027774198
      FLGFEQIFK                         0.177559156  0.585711823  0.408152666
      GANIASFVMVADAMLDQGDVF             0.458118330  1.033380564  0.575262233
      GCIISETGITSEQIHDIASAK             0.190649400  0.476203799  0.285554400
      GGLCVDLK                                   NA           NA           NA
      ICYAFMR                           0.379463866  0.186332831 -0.193131035
      NSWEGVLTGK                        0.250479727  0.284333199  0.033853472
      SLEEIVDEYSTFSESK                  0.102797277  0.005088195 -0.097709082
      VLPIVSVPER                       -0.011583129 -0.305633034 -0.294049905
      VTISGSGNVAQYAALK                           NA           NA           NA
      VTWENDNGEQEVAQGYR                -0.283890118 -0.061062035  0.222828083
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR              0.124414991  0.310261617  0.185846625
      ALVAQGVK                          0.195137769           NA           NA
      FIAEGSNMGSTPEAIAVFETAR            0.122738027  0.150012895  0.027274868
      GANIASFIK                         0.284219274  0.032694966 -0.251524307
      GCIISETGITSEQVADISSAK             0.008119322           NA           NA
      HIGQDTDVPAGDIGVGGR                0.310115657  0.432189004  0.122073348
      IMINCFNECIDYAK                    0.063851598  0.111932968  0.048081370
      ITWTSER                           0.053313078           NA           NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     0.253120569           NA           NA
      SLEQIVNEYSTFSENK                 -0.119046785 -0.642492835 -0.523446050
      STATGPSEAVWYGPPK                 -0.008255353 -0.139884005 -0.131628653
      VDIALPCATQNEVSGEEAK               0.137418487 -0.547688428 -0.685106915
      VIELGGTVVSLSDSK                   0.134703848  0.325631992  0.190928144
      VQYIAGARPWTHVQK                            NA           NA           NA
      VTWENDKGEQEVAQGYR                 0.230567742  0.711498870  0.480931128
      AAGLTAAYAR                        0.241960351  0.028851315 -0.213109037
      APEAEQVLSAAATFPIAQPATDVEAR       -0.646411615 -2.870632628 -2.224221013
      AVQDNGESAFR                       0.087845128  0.343001218  0.255156090
      DGKAPEAEQVLSAAATFPIAQPATDVEAR    -0.099604997  0.051303168  0.150908165
      GFTLAEVK                          0.238925557  0.315601206  0.076675649
      IAPRPLDLLRPVVR                    0.901416995 -0.103267139 -1.004684134
      IIVFPR                           -0.096228321 -1.696671253 -1.600442932
      NQEIFDANVQR                       0.010504533 -0.147953524 -0.158458057
      TIGIAVDHR                        -0.173148784  0.226234219  0.399383003
      VHFDQAGK                          0.297945240  0.698011315  0.400066075
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                         0.031061042  0.144856716  0.113795674
      ANGTTVLVGMPAGAK                            NA  0.070741522           NA
      ATDGGAHGVINVSVSEAAIEASTR          0.074050226  0.226375297  0.152325072
      CCSDVFNQVVK                       0.103324574  0.469767710  0.366443137
      DIVGAVLK                                   NA           NA  0.147652336
      EALDFFAR                          0.117470999  0.132154022  0.014683024
      EKDIVGAVLK                        0.146917628  0.388119775  0.241202147
      GVIFYESHGK                        0.164939215  1.216259594  1.051320379
      IGDYAGIK                         -0.072100905           NA           NA
      LPLVGGHEGAGVVVGMGENVK             0.140544229 -0.039030600 -0.179574829
      SANLMAGHWVAISGAAGGLGSLAVQYAK      0.795754496  1.273220618  0.477466122
      SIGGEVFIDFTK                      0.134709136  0.061943766 -0.072765370
      SIPETQK                          -0.376425426  0.529656517  0.906081943
      SISIVGSYVGNR                      0.064738901  0.030475636 -0.034263264
      VLGIDGGEGK                       -1.039854450 -0.991973071  0.047881379
      VLGIDGGEGKEELFR                   0.242450665  0.628276258  0.385825593
      VVGLSTLPEIYEK                     0.076880090  0.022405713 -0.054474377
      YSGVCHTDLHAWHGDWPLPVK            -0.260950377 -0.375622214 -0.114671838
      ANGTVVLVGLPAGAK                   0.277295218  0.422712134  0.145416916
      CSSDVFNHVVK                                NA           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                    -0.160224404 -0.230783404 -0.070559000
      DIPVPEPKPNEILINVK                          NA -0.210798301           NA
      EALDFFSR                          0.267374217  0.056662928 -0.210711289
      GVIFYENK                                   NA           NA  0.360868907
      IQQGTDLAEVAPILCAGVTVYK            0.080160685  0.011710403 -0.068450282
      IVGLSELPK                         0.181997462 -0.166428651 -0.348426113
      NMVSDIQEATK                      -0.009937102  0.312033525  0.321970626
      VLGIDAGEEK                                 NA           NA           NA

---

    Code
      SummarizedExperiment::rowData(D1)
    Output
      DataFrame with 87 rows and 1 column
                                           Sequence
                                        <character>
      AAADALSDLEIK                     AAADALSDLEIK
      AAADALSDLEIKDSK               AAADALSDLEIKDSK
      AEWALR                                 AEWALR
      DEGLHTDFACLLFAHLK           DEGLHTDFACLLFAHLK
      DIHDWNNR                             DIHDWNNR
      ...                                       ...
      GVIFYENK                             GVIFYENK
      IQQGTDLAEVAPILCAGVTVYK IQQGTDLAEVAPILCAGVTVYK
      IVGLSELPK                           IVGLSELPK
      NMVSDIQEATK                       NMVSDIQEATK
      VLGIDAGEEK                         VLGIDAGEEK

---

    Code
      SummarizedExperiment::colData(D1)
    Output
      DataFrame with 36 rows and 1 column
                     comparison
                    <character>
      logRatio_1_2 logRatio_1_2
      logRatio_1_3 logRatio_1_3
      logRatio_1_4 logRatio_1_4
      logRatio_1_5 logRatio_1_5
      logRatio_1_6 logRatio_1_6
      ...                   ...
      logRatio_6_8 logRatio_6_8
      logRatio_6_9 logRatio_6_9
      logRatio_7_8 logRatio_7_8
      logRatio_7_9 logRatio_7_9
      logRatio_8_9 logRatio_8_9

# normalize peptide data

    Code
      D_norm_loess
    Output
      class: SummarizedExperiment 
      dim: 87 27 
      metadata(0):
      assays(1): intensities_norm
      rownames(87): AAADALSDLEIK AAADALSDLEIKDSK ... NMVSDIQEATK VLGIDAGEEK
      rowData names(1): Sequence
      colnames(27): 12500amol_R1 12500amol_R2 ... 50amol_R2 50amol_R3
      colData names(1): sample

---

    Code
      SummarizedExperiment::rowData(D_norm_loess)
    Output
      DataFrame with 87 rows and 1 column
                                           Sequence
                                        <character>
      AAADALSDLEIK                     AAADALSDLEIK
      AAADALSDLEIKDSK               AAADALSDLEIKDSK
      AEWALR                                 AEWALR
      DEGLHTDFACLLFAHLK           DEGLHTDFACLLFAHLK
      DIHDWNNR                             DIHDWNNR
      ...                                       ...
      GVIFYENK                             GVIFYENK
      IQQGTDLAEVAPILCAGVTVYK IQQGTDLAEVAPILCAGVTVYK
      IVGLSELPK                           IVGLSELPK
      NMVSDIQEATK                       NMVSDIQEATK
      VLGIDAGEEK                         VLGIDAGEEK

---

    Code
      SummarizedExperiment::colData(D_norm_loess)
    Output
      DataFrame with 27 rows and 1 column
                         sample
                    <character>
      12500amol_R1 12500amol_R1
      12500amol_R2 12500amol_R2
      12500amol_R3 12500amol_R3
      125amol_R1     125amol_R1
      125amol_R2     125amol_R2
      ...                   ...
      500amol_R2     500amol_R2
      500amol_R3     500amol_R3
      50amol_R1       50amol_R1
      50amol_R2       50amol_R2
      50amol_R3       50amol_R3

---

    Code
      D_norm_lts
    Output
      class: SummarizedExperiment 
      dim: 87 27 
      metadata(0):
      assays(1): intensities_norm
      rownames(87): AAADALSDLEIK AAADALSDLEIKDSK ... NMVSDIQEATK VLGIDAGEEK
      rowData names(1): Sequence
      colnames(27): 12500amol_R1 12500amol_R2 ... 50amol_R2 50amol_R3
      colData names(1): sample

---

    Code
      SummarizedExperiment::rowData(D_norm_lts)
    Output
      DataFrame with 87 rows and 1 column
                                           Sequence
                                        <character>
      AAADALSDLEIK                     AAADALSDLEIK
      AAADALSDLEIKDSK               AAADALSDLEIKDSK
      AEWALR                                 AEWALR
      DEGLHTDFACLLFAHLK           DEGLHTDFACLLFAHLK
      DIHDWNNR                             DIHDWNNR
      ...                                       ...
      GVIFYENK                             GVIFYENK
      IQQGTDLAEVAPILCAGVTVYK IQQGTDLAEVAPILCAGVTVYK
      IVGLSELPK                           IVGLSELPK
      NMVSDIQEATK                       NMVSDIQEATK
      VLGIDAGEEK                         VLGIDAGEEK

---

    Code
      SummarizedExperiment::colData(D_norm_lts)
    Output
      DataFrame with 27 rows and 1 column
                         sample
                    <character>
      12500amol_R1 12500amol_R1
      12500amol_R2 12500amol_R2
      12500amol_R3 12500amol_R3
      125amol_R1     125amol_R1
      125amol_R2     125amol_R2
      ...                   ...
      500amol_R2     500amol_R2
      500amol_R3     500amol_R3
      50amol_R1       50amol_R1
      50amol_R2       50amol_R2
      50amol_R3       50amol_R3

