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
      colnames(27): 12500amol_R1 12500amol_R2 ... 00050amol_R2 00050amol_R3
      colData names(1): sample

---

    Code
      SummarizedExperiment::assays(D_norm_loess)$intensities_norm
    Output
                                       12500amol_R1 12500amol_R2 12500amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                       7416635      7422671    5204482.7
      AEWALR                                4739243      4580220    5141394.5
      DEGLHTDFACLLFAHLK                     4198307      3531891     509283.8
      DIHDWNNR                              1721905      1847182           NA
      ELETLREENR                                 NA           NA           NA
      ESEFLFNAIHTIPEIGEK                   32733791     34682498   32656039.7
      GMMPGLTFSNELICR                       8472384      7023397    8183525.5
      IVTEAVEIEQR                          18027819     18304574   18291701.3
      LLVAFGNK                             11162721     12168501   11772384.4
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                           15950962     15665825   16614709.5
      TNFFEK                                6590387      6952817    5082097.9
      TVLFPIK                              13380781     12866008   12659015.5
      VENPFDFMENISLAGK                     16266913     17401460   18821562.6
      WIQDADALFGER                         17111374     18051320   17738078.9
      YFLDALPVALLGMNADLMNQYVEFVADR         25315727     28774368   27379408.8
      AANLGGVAVSGLEMAQNSQK                  8881773      3392067    9925860.0
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           4103955           NA    4875044.1
      FLGFEQIFK                            40213371     39792867   40460533.0
      GANIASFVMVADAMLDQGDVF                57354293     55850765   57667510.7
      GCIISETGITSEQIHDIASAK                 7997314      7958432    8396915.0
      GGLCVDLK                                   NA      9835545           NA
      ICYAFMR                               8974924      9095161    9623611.5
      NSWEGVLTGK                           11753257     11406663   14079548.0
      SLEEIVDEYSTFSESK                      8271436      8109435   10243783.4
      VLPIVSVPER                           39520346     36196851   38259319.9
      VTISGSGNVAQYAALK                      4286745       138829    5815654.2
      VTWENDNGEQEVAQGYR                     3257346           NA    3429338.1
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 28258426     30279532   26757062.8
      ALVAQGVK                             10655288     12590606   12249998.1
      FIAEGSNMGSTPEAIAVFETAR               23083046     39878816   42817680.9
      GANIASFIK                            18509673     14960322   16298379.2
      GCIISETGITSEQVADISSAK                 3058336      3665458    3984126.2
      HIGQDTDVPAGDIGVGGR                   27399781     31190569   33032892.4
      IMINCFNECIDYAK                       13660929     13203073   12282716.4
      ITWTSER                              19549200     16703648   16114885.9
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        11048170     88779743   93615828.1
      SLEQIVNEYSTFSENK                     54072330     68534940   56621909.4
      STATGPSEAVWYGPPK                     48825629     39353924   39344240.4
      VDIALPCATQNEVSGEEAK                  48229476     44995617   49809682.7
      VIELGGTVVSLSDSK                      13225663     14041577   15683809.4
      VQYIAGARPWTHVQK                            NA      6506194    4237986.9
      VTWENDKGEQEVAQGYR                    33558987     31950610   37238452.8
      AAGLTAAYAR                           64343524     60923547   52881653.6
      APEAEQVLSAAATFPIAQPATDVEAR           28126195     27725671   27101281.5
      AVQDNGESAFR                          10217764     12654514   12273544.2
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        26692945     27041012   27042296.2
      GFTLAEVK                             27842784     28411089   28092123.5
      IAPRPLDLLRPVVR                       16754455     16066693   18270553.5
      IIVFPR                               49587938     46654027   48597264.1
      NQEIFDANVQR                          88021507     89785261   86509500.3
      TIGIAVDHR                            23866895     25193834   22752082.5
      VHFDQAGK                              8972222      9167903   11681731.7
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           174221810    186562319  187127305.0
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR            113803427    114339326  122749073.7
      CCSDVFNQVVK                          50275679     49425370   47688844.2
      DIVGAVLK                                   NA     51792479           NA
      EALDFFAR                            202019649    201811851  197661647.2
      EKDIVGAVLK                           37814136     36777170   40438256.0
      GVIFYESHGK                           91701122     97973997  105187602.7
      IGDYAGIK                            148363681    135914103  144477417.6
      LPLVGGHEGAGVVVGMGENVK               125376781    127586047  149353735.0
      SANLMAGHWVAISGAAGGLGSLAVQYAK         27392285     28431706   22940451.0
      SIGGEVFIDFTK                         72072672     66471713   57692734.2
      SIPETQK                              28908360     25185480   25328950.7
      SISIVGSYVGNR                        281866423    303709670  278947607.5
      VLGIDGGEGK                           20741217     19586280           NA
      VLGIDGGEGKEELFR                     318688506    332672639  340725647.9
      VVGLSTLPEIYEK                       597976972    567406017  570717500.7
      YSGVCHTDLHAWHGDWPLPVK               102581139     96101975  106810704.7
      ANGTVVLVGLPAGAK                       4736906      4366592    4477926.5
      CSSDVFNHVVK                                NA           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        38796173     35176671   36816170.0
      DIPVPEPKPNEILINVK                    17840696     20001776   16407888.8
      EALDFFSR                              7006424      5653179    6224268.1
      GVIFYENK                              5134860      5125689    5145040.1
      IQQGTDLAEVAPILCAGVTVYK                7723681      7742356    5856528.0
      IVGLSELPK                            16190476     15613610   17812692.3
      NMVSDIQEATK                           5366231      6590470    5653659.4
      VLGIDAGEEK                                 NA           NA           NA
                                       00125amol_R1 00125amol_R2 00125amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                      12542134   12748316.2   13056520.7
      AEWALR                                5873341    5540559.2    5171283.9
      DEGLHTDFACLLFAHLK                     2908910    4561543.2    4044932.8
      DIHDWNNR                              2015931           NA    2466823.1
      ELETLREENR                            5027404    4577877.3    6463092.4
      ESEFLFNAIHTIPEIGEK                   32853023   30337885.6   37945923.1
      GMMPGLTFSNELICR                       6428364    7248572.3    7801093.4
      IVTEAVEIEQR                          16794056   18063588.0   17852505.2
      LLVAFGNK                              6342165    7151192.7    7531635.7
      LLVAFGNKK                             1614911    1849261.2    1759615.8
      NKPDPAIVEK                           18535825   17715883.0   16417637.0
      TNFFEK                                5620719    4901469.8    4851619.0
      TVLFPIK                              11269724   10042297.9   10166416.1
      VENPFDFMENISLAGK                     18012043   16587488.6   14723486.8
      WIQDADALFGER                         16010634   15712144.5   16407135.3
      YFLDALPVALLGMNADLMNQYVEFVADR         24271214   21174281.0   20173803.8
      AANLGGVAVSGLEMAQNSQK                  9672333    9089502.5    9581929.1
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                           23627890   24219535.7   22615842.2
      FHPSVNLSILK                           4755650    4584150.0    4978990.4
      FLGFEQIFK                            49881220   44323364.0   45340457.0
      GANIASFVMVADAMLDQGDVF                82806757   82922348.8   77655061.2
      GCIISETGITSEQIHDIASAK                 8134665    7280455.9    8517477.3
      GGLCVDLK                             10665277           NA           NA
      ICYAFMR                              11427403   12857157.8   11645650.8
      NSWEGVLTGK                           12755482   17012304.9           NA
      SLEEIVDEYSTFSESK                      9098973    9445304.0    8879907.8
      VLPIVSVPER                           33322066   37324809.5   31087628.6
      VTISGSGNVAQYAALK                      6604325     103644.1           NA
      VTWENDNGEQEVAQGYR                          NA    2770148.7    2294475.1
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA    1794945.1
      AANLGGVAVSGLEMAQNSQR                  9312168   33211567.5   32652793.6
      ALVAQGVK                             14086194   13552945.9   13463427.6
      FIAEGSNMGSTPEAIAVFETAR               40560858   40489619.7   41542548.0
      GANIASFIK                            17916925   17932457.3   18697197.5
      GCIISETGITSEQVADISSAK                 2843776    2901471.9    2664560.8
      HIGQDTDVPAGDIGVGGR                   34078471   33868688.8   31820530.3
      IMINCFNECIDYAK                       12212156   12406984.4   12043416.6
      ITWTSER                              16821824   16336302.7     512132.4
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR         9804523  111186070.0   97218030.0
      SLEQIVNEYSTFSENK                     47249070   50363173.0   43543560.3
      STATGPSEAVWYGPPK                     38642636   42040134.9   45069896.7
      VDIALPCATQNEVSGEEAK                  47839363   51232215.6   48089298.8
      VIELGGTVVSLSDSK                      12747983   12625958.9   12732975.9
      VQYIAGARPWTHVQK                       3600738           NA    3689920.2
      VTWENDKGEQEVAQGYR                    42941014   39636484.8   42689407.5
      AAGLTAAYAR                           60589124   59640461.6   60141599.2
      APEAEQVLSAAATFPIAQPATDVEAR           18958638   23282726.0   20270461.4
      AVQDNGESAFR                          14135415   13047834.0   14247353.0
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        39055680   40605770.8   39335076.7
      GFTLAEVK                             31147043   31951191.1   30538197.3
      IAPRPLDLLRPVVR                       22326377   27021048.6   26667643.6
      IIVFPR                               40532907   42152662.6    3100005.1
      NQEIFDANVQR                          82296876   79138526.6   78484591.8
      TIGIAVDHR                            19131324   18354980.2   16612327.7
      VHFDQAGK                                   NA           NA   13324735.6
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           169563460  175449912.5  169019554.9
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             83960030   62508067.2   79746698.6
      CCSDVFNQVVK                          57673050   58238322.8   52549448.5
      DIVGAVLK                                   NA   56637988.6   49473027.4
      EALDFFAR                            204858230  203571735.8  202092626.7
      EKDIVGAVLK                           31313290   32190768.9   33096630.2
      GVIFYESHGK                          119554028  115411136.6  107276909.2
      IGDYAGIK                            138538615  133035794.9  138068558.7
      LPLVGGHEGAGVVVGMGENVK               115244883   83319484.3  147538872.9
      SANLMAGHWVAISGAAGGLGSLAVQYAK         21905762   27489334.2   31170430.0
      SIGGEVFIDFTK                         81682488   61177485.9   73589639.8
      SIPETQK                              25975869   44071024.8   45222355.1
      SISIVGSYVGNR                        263693777  281156543.2  280807539.4
      VLGIDGGEGK                            9534285   11378709.1   10427696.7
      VLGIDGGEGKEELFR                     425342325  392046557.3  337211765.9
      VVGLSTLPEIYEK                       563522051  558677990.9  567780800.9
      YSGVCHTDLHAWHGDWPLPVK                93102737  103309409.6  115080725.6
      ANGTVVLVGLPAGAK                       4528679    5662906.9    4879655.6
      CSSDVFNHVVK                           3955710    4696088.4    4096360.6
      DIPVPKPKPNELLINVK                          NA    1174684.6           NA
      VVGLSSLPEIYEK                        35533162   34835005.6   33966899.6
      DIPVPEPKPNEILINVK                    17329289   20697881.6   19867095.1
      EALDFFSR                              6566787    6993565.3    6306427.6
      GVIFYENK                                   NA    5739416.2    5542046.0
      IQQGTDLAEVAPILCAGVTVYK                7617729    7265801.8    7956923.2
      IVGLSELPK                            17847981   18546307.2   20189598.8
      NMVSDIQEATK                           3538173    4143944.7           NA
      VLGIDAGEEK                                 NA           NA           NA
                                       25000amol_R1 25000amol_R2 25000amol_R3
      AAADALSDLEIK                               NA           NA     30007758
      AAADALSDLEIKDSK                       7426422    9197437.9      5235166
      AEWALR                                4773440    5302565.9      5407099
      DEGLHTDFACLLFAHLK                     5065245    3281418.6      5832729
      DIHDWNNR                              2594907     711205.6      2238423
      ELETLREENR                                 NA    1052404.4           NA
      ESEFLFNAIHTIPEIGEK                   33193356   27787419.3     34066748
      GMMPGLTFSNELICR                       7904101    8153208.9      7074071
      IVTEAVEIEQR                          18710476   19693617.5     21844516
      LLVAFGNK                             10789786   13959025.3     11369297
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                           15012747   15057996.2     18796406
      TNFFEK                                5680397    6352153.2      7291564
      TVLFPIK                              12017867   13530087.9     12562054
      VENPFDFMENISLAGK                     17501188   15200427.7     14104791
      WIQDADALFGER                         18471130   16779078.0     16430272
      YFLDALPVALLGMNADLMNQYVEFVADR         30563930   21129679.3     21777058
      AANLGGVAVSGLEMAQNSQK                  9050105    9730415.1      8107972
      DAVWFGPPK                           178599386           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           5848132    6940842.6      6530260
      FLGFEQIFK                            38003192   37919534.8     33461565
      GANIASFVMVADAMLDQGDVF                55107486   43283425.0     44748433
      GCIISETGITSEQIHDIASAK                 7627461    6446319.9      8299955
      GGLCVDLK                             10038280   10778817.8     11163027
      ICYAFMR                               8507549    6245417.9      8551858
      NSWEGVLTGK                           13681165   13838803.9     12493499
      SLEEIVDEYSTFSESK                      9954196    6805743.9      8440860
      VLPIVSVPER                           40312059   49996703.2     39852430
      VTISGSGNVAQYAALK                      1901774    5943212.1      1096516
      VTWENDNGEQEVAQGYR                     2515030           NA      2014214
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 31051078   29853067.3     26691774
      ALVAQGVK                             12249580   16102486.3     16163001
      FIAEGSNMGSTPEAIAVFETAR               39274882   42043886.9     36057121
      GANIASFIK                            17326865   16084564.1           NA
      GCIISETGITSEQVADISSAK                 3010007    3180561.9      3434536
      HIGQDTDVPAGDIGVGGR                   30139290   33135486.3     34592672
      IMINCFNECIDYAK                       11268616   10934682.1     10579596
      ITWTSER                              17233894   20112809.9     19450227
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR       101512428   70738820.6     84653312
      SLEQIVNEYSTFSENK                     56553823   50830814.0     55808537
      STATGPSEAVWYGPPK                     42686060   50202589.7     40237527
      VDIALPCATQNEVSGEEAK                  50987256   53788881.4     48677212
      VIELGGTVVSLSDSK                      13644298   14719650.6     13825702
      VQYIAGARPWTHVQK                       3570062    6313093.9      2809371
      VTWENDKGEQEVAQGYR                    34454808   30047812.5     30641659
      AAGLTAAYAR                           59305094   68204319.9     68398379
      APEAEQVLSAAATFPIAQPATDVEAR           30719427   27609674.2     30804089
      AVQDNGESAFR                          11551026   12763624.0     13559097
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        18993936   22212858.2     22539664
      GFTLAEVK                             30002402   27425480.8     26779955
      IAPRPLDLLRPVVR                       16714495   13054293.6     17747871
      IIVFPR                               51811392   58366768.3     50366002
      NQEIFDANVQR                          83986018   91819323.6     90192053
      TIGIAVDHR                            24773682   23835943.0     27115740
      VHFDQAGK                             11652336   10168391.1     11375885
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           190336933  198574368.9    188960926
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             86600534   90842065.5    115142127
      CCSDVFNQVVK                          47393922   51163574.1     47505751
      DIVGAVLK                             49538299   57652599.1     50254912
      EALDFFAR                            204340463  225439818.6    191509885
      EKDIVGAVLK                           38606091   31155015.9     30687542
      GVIFYESHGK                           99001930   95710930.8    100353621
      IGDYAGIK                            152015207  180150042.8    177891441
      LPLVGGHEGAGVVVGMGENVK               130009911   99721363.6    121335243
      SANLMAGHWVAISGAAGGLGSLAVQYAK         28081118   41871149.2     31272385
      SIGGEVFIDFTK                         94187179   58149712.2     50080502
      SIPETQK                              23700810   19966557.6     13597713
      SISIVGSYVGNR                        295391281  305309426.9    300044196
      VLGIDGGEGK                           19896071   26916786.9     27706308
      VLGIDGGEGKEELFR                     330366360  289103294.0    302426716
      VVGLSTLPEIYEK                       551215027  593571143.0    567316893
      YSGVCHTDLHAWHGDWPLPVK               117551726   90737588.2    111419963
      ANGTVVLVGLPAGAK                       4369083    5037733.4      4712681
      CSSDVFNHVVK                           4268312           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        40265595   41258366.9     37419933
      DIPVPEPKPNEILINVK                    21017884   21923942.6     21244884
      EALDFFSR                              6305159    6699293.3      4995221
      GVIFYENK                              4585439    4197136.5      6020083
      IQQGTDLAEVAPILCAGVTVYK                8037603    6879064.2      6018069
      IVGLSELPK                            18481767   16444396.2     15768297
      NMVSDIQEATK                           5713499    6490342.2      6125405
      VLGIDAGEEK                                 NA           NA           NA
                                       02500amol_R1 02500amol_R2 02500amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                     5595921.4    7304945.4      8271873
      AEWALR                              5415551.0    4585275.8      4617641
      DEGLHTDFACLLFAHLK                   4768539.0    4990409.6      3619076
      DIHDWNNR                             585271.2           NA           NA
      ELETLREENR                                 NA           NA      4393575
      ESEFLFNAIHTIPEIGEK                 34182940.1   26125178.6     35222206
      GMMPGLTFSNELICR                     8147712.1    9354975.0      8960250
      IVTEAVEIEQR                        17267374.5   17738445.1     18060722
      LLVAFGNK                           10508926.9   10237990.8      9197303
      LLVAFGNKK                            709934.6           NA           NA
      NKPDPAIVEK                         14143215.8   14219062.6     17624177
      TNFFEK                              6460702.8    5614478.2      5865591
      TVLFPIK                             9965079.5   11222433.1     11864199
      VENPFDFMENISLAGK                   16396797.8   15349285.9     15770832
      WIQDADALFGER                       18833245.5   21082462.1     18654488
      YFLDALPVALLGMNADLMNQYVEFVADR       24372178.0   23197181.4     24352055
      AANLGGVAVSGLEMAQNSQK               10490770.7   10176421.3      8081602
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                         5278015.6    4647559.9      7440398
      FLGFEQIFK                          40431797.2   43213593.9     39087378
      GANIASFVMVADAMLDQGDVF              52754092.0   58689985.4     57699983
      GCIISETGITSEQIHDIASAK               5926437.0    7847804.5      8237986
      GGLCVDLK                            9530173.8    9183527.3     10902355
      ICYAFMR                             8255233.6    8583843.8      9412257
      NSWEGVLTGK                         12440598.2   13039294.0     11450575
      SLEEIVDEYSTFSESK                    8167129.2   10654697.0      7994925
      VLPIVSVPER                         41327846.7   40053853.2     42560449
      VTISGSGNVAQYAALK                    6348398.1     183734.1           NA
      VTWENDNGEQEVAQGYR                   3214988.8    2899980.2      3027215
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR               31389393.2   33629110.8     32491830
      ALVAQGVK                                   NA   12119030.7     13000896
      FIAEGSNMGSTPEAIAVFETAR             36244891.0   37714643.0     35885608
      GANIASFIK                          17560953.5   17145744.4     17998722
      GCIISETGITSEQVADISSAK               3576505.6    3175210.2      3414993
      HIGQDTDVPAGDIGVGGR                 33945545.5   28525212.3     31844467
      IMINCFNECIDYAK                     12691390.1   13342830.2     13221805
      ITWTSER                                    NA   16483821.3     16273348
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR      92863232.6   99962652.5     85105749
      SLEQIVNEYSTFSENK                   53892682.9   54143103.9     58870074
      STATGPSEAVWYGPPK                   49661750.0   50743385.0     44661768
      VDIALPCATQNEVSGEEAK                54518663.4   58040234.0     54199520
      VIELGGTVVSLSDSK                    13502507.5   13816132.6     13275486
      VQYIAGARPWTHVQK                     7946852.8           NA      6059156
      VTWENDKGEQEVAQGYR                  34110502.5   32486380.1     33585567
      AAGLTAAYAR                         56954601.9   57501302.3     56560919
      APEAEQVLSAAATFPIAQPATDVEAR         25884481.2   27266448.8     28354477
      AVQDNGESAFR                        12940670.2   12659915.9     13179089
      DGKAPEAEQVLSAAATFPIAQPATDVEAR      27633187.5   22281422.5     20863834
      GFTLAEVK                           29527347.5   31537179.5     29774904
      IAPRPLDLLRPVVR                      8569716.5   19158263.6      9922301
      IIVFPR                             59040933.0   56681087.1     51072227
      NQEIFDANVQR                        91240133.0   91421926.1     88874121
      TIGIAVDHR                          21794063.5   22576316.2     23003411
      VHFDQAGK                            8931563.3           NA     10515199
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                         170439308.1  156784515.3    160751087
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR           76204965.9   97289521.2     91899865
      CCSDVFNQVVK                        53568931.3   53126544.8     53491110
      DIVGAVLK                           48747466.7           NA     49328956
      EALDFFAR                          216714599.0  212983592.6    206743632
      EKDIVGAVLK                         29432759.3   30770369.1     36041371
      GVIFYESHGK                         91540027.8   91461313.5     94153466
      IGDYAGIK                          149374178.1  136928351.8    147930090
      LPLVGGHEGAGVVVGMGENVK             109539068.8  122284859.7    158355967
      SANLMAGHWVAISGAAGGLGSLAVQYAK       36875060.9   21078521.2     28388010
      SIGGEVFIDFTK                       54301505.8   74023463.4     73768407
      SIPETQK                            39192219.3   21023741.0     23565117
      SISIVGSYVGNR                      297295059.2  285883550.2    272141962
      VLGIDGGEGK                         18885089.3   18871600.2     20616709
      VLGIDGGEGKEELFR                   319004938.9  336174591.8    354519196
      VVGLSTLPEIYEK                     619743561.9  577591524.8    548306769
      YSGVCHTDLHAWHGDWPLPVK             157958195.9  160081306.9    146633703
      ANGTVVLVGLPAGAK                     4774774.0    3874026.4      4891201
      CSSDVFNHVVK                                NA           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                      36342577.4   38410389.8     34364001
      DIPVPEPKPNEILINVK                  21340108.1   22893162.2     21609988
      EALDFFSR                            6914459.3    6235909.3      8525426
      GVIFYENK                            4784954.9    5065378.2      5453076
      IQQGTDLAEVAPILCAGVTVYK              6450726.9    8199292.0      6737493
      IVGLSELPK                          18687714.7   18910780.4     19273889
      NMVSDIQEATK                         6041523.4    6643690.1      5890578
      VLGIDAGEEK                                 NA           NA           NA
                                       00250amol_R1 00250amol_R2 00250amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                    10432806.0     13449987   12305697.7
      AEWALR                              4992311.8      5973545    4944225.4
      DEGLHTDFACLLFAHLK                          NA           NA           NA
      DIHDWNNR                                   NA      1784085    1989942.3
      ELETLREENR                          3995045.4      3908644    3969456.0
      ESEFLFNAIHTIPEIGEK                 27604726.8     24326285   35914279.8
      GMMPGLTFSNELICR                     7487400.0           NA    7587894.8
      IVTEAVEIEQR                        18064901.1     16707881   18751362.6
      LLVAFGNK                            6908938.5      6293379    7254632.3
      LLVAFGNKK                                  NA      1699510           NA
      NKPDPAIVEK                         16426133.5     19258645   18373233.2
      TNFFEK                              6535226.7      5580499    5342732.1
      TVLFPIK                             9303627.3      9618992   11579265.7
      VENPFDFMENISLAGK                   17512067.9     17502116   17872966.6
      WIQDADALFGER                       14461339.6     14556023   16818578.8
      YFLDALPVALLGMNADLMNQYVEFVADR       23849756.7     22783680   23454879.2
      AANLGGVAVSGLEMAQNSQK               10631235.0     10185417    5929991.8
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA     20403251   20988333.1
      FHPSVNLSILK                         4720430.6      6412563    5273298.5
      FLGFEQIFK                          47019736.6     45764537   44180953.5
      GANIASFVMVADAMLDQGDVF              76947303.2     86545477   79085293.6
      GCIISETGITSEQIHDIASAK               8560236.6      8729617    6946353.8
      GGLCVDLK                           11369315.4     11276818   10656253.9
      ICYAFMR                            10277775.0     10629209   12138396.9
      NSWEGVLTGK                         12212103.1     13669251   14091650.2
      SLEEIVDEYSTFSESK                    8621836.4      8114059    8135521.2
      VLPIVSVPER                         33329281.1     33504649   39131442.3
      VTISGSGNVAQYAALK                    3104425.9           NA     146341.6
      VTWENDNGEQEVAQGYR                   2801033.4      2622178    3050225.3
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR               36431866.7     34044441   38060786.4
      ALVAQGVK                           13791725.6     12756835   12104892.4
      FIAEGSNMGSTPEAIAVFETAR             42432836.1     27495977   42053267.3
      GANIASFIK                          17612768.2     17725728   18707577.4
      GCIISETGITSEQVADISSAK               3993055.1      3569734    3378317.7
      HIGQDTDVPAGDIGVGGR                 33021732.7     36498814   30742386.8
      IMINCFNECIDYAK                     14410944.4     13056425   13080203.7
      ITWTSER                            17379144.5     16604256   16383670.7
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     107241136.7    101560420  100124642.2
      SLEQIVNEYSTFSENK                   54617831.1     54585265   60289714.9
      STATGPSEAVWYGPPK                   46680169.6     45631617   38012429.1
      VDIALPCATQNEVSGEEAK                49060041.3     43755088   50772575.1
      VIELGGTVVSLSDSK                    13619572.3     12941267   13044742.3
      VQYIAGARPWTHVQK                    10381297.8      9863001    4480124.5
      VTWENDKGEQEVAQGYR                  36355856.8     42159340   38569705.4
      AAGLTAAYAR                         63970386.8     59045494   59896802.9
      APEAEQVLSAAATFPIAQPATDVEAR         20009397.1     17922721   14500108.2
      AVQDNGESAFR                          274611.5     13563451   11655732.7
      DGKAPEAEQVLSAAATFPIAQPATDVEAR      27156709.1     26512344   28095155.5
      GFTLAEVK                           30428215.0     33973687   30238177.0
      IAPRPLDLLRPVVR                     12899141.8     23735298   23649149.3
      IIVFPR                             43159966.6     48003990   42803039.6
      NQEIFDANVQR                        93541636.4     85582423   80174676.8
      TIGIAVDHR                          20962741.5     22560720   20002487.8
      VHFDQAGK                           10766230.8     10677440   10312928.0
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                         173667438.7    181977833  166237723.1
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR           80449913.2     80408883   67675191.6
      CCSDVFNQVVK                        52290145.6     58027488   56283831.6
      DIVGAVLK                           54119455.0     56428952           NA
      EALDFFAR                          203124470.9    198394114  214770566.1
      EKDIVGAVLK                         33747669.4     36613552   36423422.8
      GVIFYESHGK                        103354891.0     95766230  103263799.6
      IGDYAGIK                          140895539.2    137937565           NA
      LPLVGGHEGAGVVVGMGENVK             124064933.3    128741667   88566915.2
      SANLMAGHWVAISGAAGGLGSLAVQYAK       26581976.2     30044570   28405627.7
      SIGGEVFIDFTK                       65183511.7     62400220   54518382.1
      SIPETQK                            27720338.2     23181190   24266755.2
      SISIVGSYVGNR                      264778952.2    277006490  276916286.5
      VLGIDGGEGK                         10002023.0     10604452   10476591.2
      VLGIDGGEGKEELFR                   345005802.9    363059894  372389220.0
      VVGLSTLPEIYEK                     595091862.0    558268160  578392399.9
      YSGVCHTDLHAWHGDWPLPVK             103108222.4    118396106  113011945.5
      ANGTVVLVGLPAGAK                     4586792.5      4539262    5219552.7
      CSSDVFNHVVK                                NA      3753272    3936969.3
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                      36671690.0     33881009   33305631.6
      DIPVPEPKPNEILINVK                    754910.9     19071471   20663071.0
      EALDFFSR                            7234107.8           NA    6844113.6
      GVIFYENK                                   NA      5466692    5916488.1
      IQQGTDLAEVAPILCAGVTVYK              9466923.4      6288807    9109788.0
      IVGLSELPK                          17726354.2     17824261   19450720.4
      NMVSDIQEATK                         5781644.8      5349892    6317593.9
      VLGIDAGEEK                                 NA           NA           NA
                                       50000amol_R1 50000amol_R2 50000amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                       3844563      2663432      9856128
      AEWALR                                4461982      4846864      4793571
      DEGLHTDFACLLFAHLK                          NA      4847097      2361907
      DIHDWNNR                              2521743      1343634           NA
      ELETLREENR                            4767176      4916800      5411949
      ESEFLFNAIHTIPEIGEK                   34687669     34300760     30387665
      GMMPGLTFSNELICR                       8050453      7406513      5773522
      IVTEAVEIEQR                          21911458     23522830     21871457
      LLVAFGNK                             14190713     15144345     13399637
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                           19029930     22130420     21508465
      TNFFEK                                4878597      5812238      5337454
      TVLFPIK                              14160097     14140119     13587674
      VENPFDFMENISLAGK                     11865567     12729190     12879355
      WIQDADALFGER                         13343425     12553696     15718482
      YFLDALPVALLGMNADLMNQYVEFVADR         22536509     23334915     23600183
      AANLGGVAVSGLEMAQNSQK                  8971557      8753020      8253881
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                          11507153      5874332      7361643
      FLGFEQIFK                            36346906     38890487     37598701
      GANIASFVMVADAMLDQGDVF                40414610     39898448     43082960
      GCIISETGITSEQIHDIASAK                 6839864      7518247      8920209
      GGLCVDLK                             10708954     12707534     10973045
      ICYAFMR                               8403346      8524347      8342594
      NSWEGVLTGK                           14325337     14138224     13260352
      SLEEIVDEYSTFSESK                      5057614      6856776      5316528
      VLPIVSVPER                           41665384     40811797     36764278
      VTISGSGNVAQYAALK                      5543059      5485772           NA
      VTWENDNGEQEVAQGYR                          NA      2663865           NA
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 28001755     22244406     23890787
      ALVAQGVK                             14286823     14133488     13720141
      FIAEGSNMGSTPEAIAVFETAR               18099736     28524224     32407822
      GANIASFIK                            17631820     18463061     18176920
      GCIISETGITSEQVADISSAK                 3387061      3998465      4552468
      HIGQDTDVPAGDIGVGGR                   39981849     43699500     35756672
      IMINCFNECIDYAK                       13908668     10777073     10446821
      ITWTSER                              14827470     16703133           NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        71916583      9738282      9500198
      SLEQIVNEYSTFSENK                     51299905     45474992     53319074
      STATGPSEAVWYGPPK                     38618555           NA     39472267
      VDIALPCATQNEVSGEEAK                  51233000     48781764     43881337
      VIELGGTVVSLSDSK                      12917519     12286762     12759434
      VQYIAGARPWTHVQK                       3825160      4344547           NA
      VTWENDKGEQEVAQGYR                    32043415     37751952     38871267
      AAGLTAAYAR                           61696195     62454344     65830130
      APEAEQVLSAAATFPIAQPATDVEAR           25296650     29340919     25964405
      AVQDNGESAFR                          12615970     12441246     14326139
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        12895952     11480780     14213356
      GFTLAEVK                             29014921     28502524     32423232
      IAPRPLDLLRPVVR                       17889149     10664324     11395310
      IIVFPR                               54107389     53866360     50012277
      NQEIFDANVQR                          79681325     84132026     87131600
      TIGIAVDHR                                  NA     23826107     26617151
      VHFDQAGK                             17520079     21831006     21167700
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           222503396    213599592    212844433
      ANGTTVLVGMPAGAK                      19338108     18105650     19231825
      ATDGGAHGVINVSVSEAAIEASTR             98064020    119794025    152149267
      CCSDVFNQVVK                          47462025     51102709     49211376
      DIVGAVLK                             53263269     52983113     54600226
      EALDFFAR                            196536271    193552033    192294131
      EKDIVGAVLK                           38813823     42118177     39277967
      GVIFYESHGK                          109790797    118810146    135567783
      IGDYAGIK                            147175599    149566572    153638830
      LPLVGGHEGAGVVVGMGENVK               136943074    142318627    119039259
      SANLMAGHWVAISGAAGGLGSLAVQYAK         18016017     18133263     35076013
      SIGGEVFIDFTK                         52829655     58533965     45370446
      SIPETQK                              11946977      9609085      8979866
      SISIVGSYVGNR                        307441900    305149356    295881240
      VLGIDGGEGK                           21011656     23812826     26372475
      VLGIDGGEGKEELFR                     334492502    334694407    317548059
      VVGLSTLPEIYEK                       529744573    474151094    526841585
      YSGVCHTDLHAWHGDWPLPVK               102967957     91637653    106519257
      ANGTVVLVGLPAGAK                       5394400      6031044      3750723
      CSSDVFNHVVK                           3961232      3840679      3776567
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        32588470     34646821     35865758
      DIPVPEPKPNEILINVK                    18934290     18933404     17770481
      EALDFFSR                              6599570      7067314           NA
      GVIFYENK                                   NA      6798971      4386911
      IQQGTDLAEVAPILCAGVTVYK                6313893      6690293      5888141
      IVGLSELPK                            19294297     18309375     19061374
      NMVSDIQEATK                           4911099      4143744      5872798
      VLGIDAGEEK                                 NA           NA           NA
                                       05000amol_R1 05000amol_R2 05000amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                            NA      2159461           NA
      AEWALR                                5338363      4831952      5118987
      DEGLHTDFACLLFAHLK                     4920523           NA           NA
      DIHDWNNR                                   NA           NA           NA
      ELETLREENR                                 NA      4392582           NA
      ESEFLFNAIHTIPEIGEK                   27644709     31189816     33679615
      GMMPGLTFSNELICR                      11147893      9664644      8754728
      IVTEAVEIEQR                          15523554     15872595     17781794
      LLVAFGNK                              9889801           NA     11068318
      LLVAFGNKK                                  NA      1452545           NA
      NKPDPAIVEK                           15072771     15671263     17344443
      TNFFEK                                5940708      7276363      6090163
      TVLFPIK                              11062400     11552145     10622215
      VENPFDFMENISLAGK                     18137691     15838824     18067055
      WIQDADALFGER                         18527088     19984554     19732356
      YFLDALPVALLGMNADLMNQYVEFVADR         21941924     21749565     23206451
      AANLGGVAVSGLEMAQNSQK                  8516938      5617373      9591444
      DAVWFGPPK                           170061486           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           5078162      4702085      4946042
      FLGFEQIFK                            39933599     41978044     38171348
      GANIASFVMVADAMLDQGDVF                53533682     49417035     54030561
      GCIISETGITSEQIHDIASAK                 6956225      7399136      6721639
      GGLCVDLK                             10474020           NA      8895072
      ICYAFMR                               9114244      9443498      8932356
      NSWEGVLTGK                           10204301     11703872     12493455
      SLEEIVDEYSTFSESK                      8387129      7952363      8656890
      VLPIVSVPER                           42074886     43362604     42104625
      VTISGSGNVAQYAALK                           NA      6783340      3084721
      VTWENDNGEQEVAQGYR                     3395577      3343718      3618167
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 34140019     29222653      8785524
      ALVAQGVK                             11161206     12536802     12253549
      FIAEGSNMGSTPEAIAVFETAR               38937270     37363624     37222472
      GANIASFIK                            18184035     15629187     18485036
      GCIISETGITSEQVADISSAK                 4002455      4141917      2713149
      HIGQDTDVPAGDIGVGGR                   31374185     29975573     27218360
      IMINCFNECIDYAK                       14231396     11198217     13269747
      ITWTSER                              16877359     20203019     17993334
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        80934512     91043625     79910865
      SLEQIVNEYSTFSENK                     61233974     56465131     54755694
      STATGPSEAVWYGPPK                     42873260     46405422     49315492
      VDIALPCATQNEVSGEEAK                  51322939     47302877     54748251
      VIELGGTVVSLSDSK                      13096916     13728547     13076098
      VQYIAGARPWTHVQK                      12775162      4203795      9483316
      VTWENDKGEQEVAQGYR                    33157637     34169642     34585667
      AAGLTAAYAR                           55616495     58283717     57007655
      APEAEQVLSAAATFPIAQPATDVEAR           26817588     25947182     27427393
      AVQDNGESAFR                          13487853     13565367     12303196
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        30268397     26463898     26940080
      GFTLAEVK                             26890702     28835419     28220785
      IAPRPLDLLRPVVR                       18343004     15974788      9436019
      IIVFPR                               52138542     51804854     54451559
      NQEIFDANVQR                          97738954     88207113     95314269
      TIGIAVDHR                            26542836     26540388     24761085
      VHFDQAGK                              8910317      8896523      9333336
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           171613329    171372927    168591121
      ANGTTVLVGMPAGAK                      21027956     17226913     20029391
      ATDGGAHGVINVSVSEAAIEASTR             84686409     89323876     99008515
      CCSDVFNQVVK                          48159036     48754121     49913196
      DIVGAVLK                                   NA     52550005           NA
      EALDFFAR                            207974213    207506131    201333900
      EKDIVGAVLK                           30103477     32577075     34348787
      GVIFYESHGK                           90571563     87310166     96447658
      IGDYAGIK                            143145509    147892333    135951888
      LPLVGGHEGAGVVVGMGENVK               129429549    130279558    133121409
      SANLMAGHWVAISGAAGGLGSLAVQYAK         23752182     13951521     22650572
      SIGGEVFIDFTK                         65370571     63571722     59807101
      SIPETQK                              32783467     31673888     17939698
      SISIVGSYVGNR                        289781038    286971073    286400077
      VLGIDGGEGK                           18442172     19724086     20391014
      VLGIDGGEGKEELFR                     322427423    321458113    319942348
      VVGLSTLPEIYEK                       594869126    576903614    584625672
      YSGVCHTDLHAWHGDWPLPVK               125943753    137018512    153144561
      ANGTVVLVGLPAGAK                       3798768      4789004      4232723
      CSSDVFNHVVK                                NA      3584572      1344908
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        38224321     39105825     37931315
      DIPVPEPKPNEILINVK                    18790770     20049663     20321776
      EALDFFSR                              6257447      6531772      6538278
      GVIFYENK                                   NA           NA           NA
      IQQGTDLAEVAPILCAGVTVYK                6959755      8676757      8003591
      IVGLSELPK                            17491069     18993025     18508100
      NMVSDIQEATK                           5020829      5722995      7080694
      VLGIDAGEEK                                 NA           NA      4417712
                                       00500amol_R1 00500amol_R2 00500amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                      12028894      3695220    3200005.2
      AEWALR                                5205516      5433118    4915607.0
      DEGLHTDFACLLFAHLK                          NA           NA           NA
      DIHDWNNR                              1672862           NA           NA
      ELETLREENR                            3929625      3543355    5832267.7
      ESEFLFNAIHTIPEIGEK                   30871978     32193819   30555955.4
      GMMPGLTFSNELICR                       9854919      8153039    7608285.3
      IVTEAVEIEQR                          18146990     16209063   16270514.6
      LLVAFGNK                              7101724      8568095    8705231.0
      LLVAFGNKK                             1975394           NA           NA
      NKPDPAIVEK                           17924395     18754663   19149295.9
      TNFFEK                                6414219      5482565    6010263.2
      TVLFPIK                               9247222      9370435   11604345.1
      VENPFDFMENISLAGK                     16985029     13998349   16022462.2
      WIQDADALFGER                         16297570     14289052   14093355.5
      YFLDALPVALLGMNADLMNQYVEFVADR         26142508     23445925   21933617.9
      AANLGGVAVSGLEMAQNSQK                 10032127      9422422    9980888.0
      DAVWFGPPK                                  NA           NA  169687136.4
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           4434023      5234748    4631292.6
      FLGFEQIFK                            41143465     45642793   42086408.7
      GANIASFVMVADAMLDQGDVF                68204761     68089289   68626914.3
      GCIISETGITSEQIHDIASAK                 6297150      8966561    8007203.4
      GGLCVDLK                             10028500      9660436           NA
      ICYAFMR                              10930233     12433614   10979626.5
      NSWEGVLTGK                           13015245     11735689   13985862.7
      SLEEIVDEYSTFSESK                      9455710      8099877    8242714.9
      VLPIVSVPER                           39244665     41774348   38616358.3
      VTISGSGNVAQYAALK                           NA           NA     119941.4
      VTWENDNGEQEVAQGYR                     2801329      2717685    2737222.2
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 32547739      8676010   33845333.8
      ALVAQGVK                             12283056     13713856   13207130.4
      FIAEGSNMGSTPEAIAVFETAR               37917240     39807594   38987168.9
      GANIASFIK                            19670134     19751042   20571399.5
      GCIISETGITSEQVADISSAK                 3609093      3574961    3449278.4
      HIGQDTDVPAGDIGVGGR                   33341483     38808743   31340597.8
      IMINCFNECIDYAK                       12595901     12912685   12883808.5
      ITWTSER                              17128633     18178710   18638650.1
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        91629833     95782375   97984175.9
      SLEQIVNEYSTFSENK                     52861338     49362534   48932114.8
      STATGPSEAVWYGPPK                     42103474     46320483   42210528.3
      VDIALPCATQNEVSGEEAK                  55732974     52979994   51502488.9
      VIELGGTVVSLSDSK                      14271017     14545714   12651904.7
      VQYIAGARPWTHVQK                       5740645      9044110           NA
      VTWENDKGEQEVAQGYR                    38204917     38036568   36722639.2
      AAGLTAAYAR                           62869265     65089915   64331925.7
      APEAEQVLSAAATFPIAQPATDVEAR           18498428     16784103   13055130.7
      AVQDNGESAFR                          12118958     13694032   14119114.8
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        22879439     26225582   24369495.1
      GFTLAEVK                             29394024     31678886   32067651.5
      IAPRPLDLLRPVVR                       28743797     22431015   26371998.2
      IIVFPR                               46621283     47183038   45898093.9
      NQEIFDANVQR                          93052507     88375570   86084113.4
      TIGIAVDHR                            21749330     20759290   22519523.0
      VHFDQAGK                             10951542     10766545   10056154.6
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           176761890    164933839  148854228.4
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             81977224     92593962   96561126.7
      CCSDVFNQVVK                          46923989     53398300   49850369.2
      DIVGAVLK                             56861242     46752039   47790014.4
      EALDFFAR                            209378149    212381373  206287246.5
      EKDIVGAVLK                           35335969     34320002   31591705.1
      GVIFYESHGK                           96395374     98142272   94937990.5
      IGDYAGIK                            126342079    123236783  132210709.0
      LPLVGGHEGAGVVVGMGENVK               142249174    118624794  146332479.7
      SANLMAGHWVAISGAAGGLGSLAVQYAK         37324850     27503169   33730117.9
      SIGGEVFIDFTK                         54292690     67382788   75544500.2
      SIPETQK                              23287009     18864956   17976132.2
      SISIVGSYVGNR                        283067664    285832977  278122911.0
      VLGIDGGEGK                            9792550      8035368    9045383.8
      VLGIDGGEGKEELFR                     338101969    364863579  366795540.9
      VVGLSTLPEIYEK                       556472856    608188915  569776010.7
      YSGVCHTDLHAWHGDWPLPVK               123262907     88592963  114740819.0
      ANGTVVLVGLPAGAK                       4485421      5763864    4816405.1
      CSSDVFNHVVK                           3937093           NA    1048304.0
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        33979710     30866768   32502324.4
      DIPVPEPKPNEILINVK                          NA     18374084   20263697.0
      EALDFFSR                              6203850      6615591    9685888.8
      GVIFYENK                              4800182      5050089    5882250.5
      IQQGTDLAEVAPILCAGVTVYK                8209602      8276967    7569377.1
      IVGLSELPK                            18390483     20517818   19885671.2
      NMVSDIQEATK                           5711214      6677965    4650333.2
      VLGIDAGEEK                                 NA           NA           NA
                                       00050amol_R1 00050amol_R2 00050amol_R3
      AAADALSDLEIK                         17699581           NA           NA
      AAADALSDLEIKDSK                       1632247           NA           NA
      AEWALR                                3923771    4559219.4      4644144
      DEGLHTDFACLLFAHLK                     8834755    6156194.6      5060039
      DIHDWNNR                                   NA    1384860.1      1549919
      ELETLREENR                            5634223    6219788.3      6371415
      ESEFLFNAIHTIPEIGEK                   31189072   33078266.2     26765153
      GMMPGLTFSNELICR                       4762460           NA           NA
      IVTEAVEIEQR                          27047742   22673883.7     21152595
      LLVAFGNK                              6611060    5859788.7      6884724
      LLVAFGNKK                             1105307    1736240.4      2506560
      NKPDPAIVEK                           20042510   16136907.2     20390030
      TNFFEK                                4432844    5607903.8      5082510
      TVLFPIK                               8881226    9085972.5     10433599
      VENPFDFMENISLAGK                     12885065   13401072.9     14447306
      WIQDADALFGER                         15198594   12981732.9     14555151
      YFLDALPVALLGMNADLMNQYVEFVADR         27880550   28481760.4     26954787
      AANLGGVAVSGLEMAQNSQK                       NA    9564492.4      7760272
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                           27382534   22975623.3     26081866
      FHPSVNLSILK                          10035129    6560135.2      7111454
      FLGFEQIFK                            52755751   56919902.2     52864872
      GANIASFVMVADAMLDQGDVF               107327902   92717886.2     85188603
      GCIISETGITSEQIHDIASAK                 8309375    9114635.8      7161393
      GGLCVDLK                             10467488    9463963.0     10639763
      ICYAFMR                               9088356    8993172.9      9300292
      NSWEGVLTGK                           11046467   12822034.2     13934322
      SLEEIVDEYSTFSESK                      8710247    6855578.3      5525524
      VLPIVSVPER                           29258969   31122618.6     32810995
      VTISGSGNVAQYAALK                      2016368    2031836.1           NA
      VTWENDNGEQEVAQGYR                     2017881    2137270.7      2627531
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK      6219605    4502008.1      4173628
      AANLGGVAVSGLEMAQNSQR                 27177911   27181563.2     28487505
      ALVAQGVK                                   NA     638395.7           NA
      FIAEGSNMGSTPEAIAVFETAR               35668277   39843800.4     38465311
      GANIASFIK                            15769398   16946333.5     17316726
      GCIISETGITSEQVADISSAK                 4961946    3928327.6           NA
      HIGQDTDVPAGDIGVGGR                   39060567   36119306.9     32768624
      IMINCFNECIDYAK                       11895945   13068642.9     12689270
      ITWTSER                                    NA   17541389.0     16772904
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        82367675           NA     97214775
      SLEQIVNEYSTFSENK                     43367618    9691483.5     44552783
      STATGPSEAVWYGPPK                     37218230   38552705.8     37459610
      VDIALPCATQNEVSGEEAK                   7250299   45149446.6     45377415
      VIELGGTVVSLSDSK                      17893480   14178898.4     13173595
      VQYIAGARPWTHVQK                            NA           NA           NA
      VTWENDKGEQEVAQGYR                    52358881   50098290.7     47732593
      AAGLTAAYAR                           44931282   57976807.2     54053359
      APEAEQVLSAAATFPIAQPATDVEAR            4699123    2674199.3      2882248
      AVQDNGESAFR                          15183081   15561116.5     13800429
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        28206052   22512124.3     29094428
      GFTLAEVK                             28379631   33978906.1     32937468
      IAPRPLDLLRPVVR                       18725988    9605888.7      9656827
      IIVFPR                                1810101    2006046.0     41592096
      NQEIFDANVQR                          66861303   80829764.0     78183930
      TIGIAVDHR                            30013206   29193696.5     24895548
      VHFDQAGK                             11824961   13428604.9     14169769
      VHFDQAGKK                             1336577    1305778.8      1274760
      ANELLINVK                           158346273  172460604.7    171151880
      ANGTTVLVGMPAGAK                      19988508   18343588.0     19086935
      ATDGGAHGVINVSVSEAAIEASTR            102875121   95175502.3     83483846
      CCSDVFNQVVK                          57910096   61956354.5     63075208
      DIVGAVLK                             48417859   56859591.1     53613020
      EALDFFAR                            196283673  204115057.7    199664277
      EKDIVGAVLK                           50407706   32260639.3     30476845
      GVIFYESHGK                          190729440  214736765.8    158126419
      IGDYAGIK                                   NA  123929496.2    119625538
      LPLVGGHEGAGVVVGMGENVK               107288585  104590990.7    126452510
      SANLMAGHWVAISGAAGGLGSLAVQYAK         50704211   41479692.7     39986739
      SIGGEVFIDFTK                         72705484   41439168.4     59704027
      SIPETQK                              28432494   28770151.4     53938312
      SISIVGSYVGNR                        265148605  261369746.4    258257299
      VLGIDGGEGK                            9773961    8249920.7      8923691
      VLGIDGGEGKEELFR                     408360866  457542691.4    467001271
      VVGLSTLPEIYEK                       526751476  533621497.9    540791821
      YSGVCHTDLHAWHGDWPLPVK                70809305  117585909.5     98142759
      ANGTVVLVGLPAGAK                       4210887    4340998.0      4545670
      CSSDVFNHVVK                           3987211    3123477.2      4281769
      DIPVPKPKPNELLINVK                     5299097    4679747.6      3518777
      VVGLSSLPEIYEK                        27620701   29100834.8     32191543
      DIPVPEPKPNEILINVK                     3755051   24601995.6     19241921
      EALDFFSR                              5760488    5421135.1      5143937
      GVIFYENK                              5071227    5726503.8      5560622
      IQQGTDLAEVAPILCAGVTVYK                5873261    6788699.9      7193493
      IVGLSELPK                            14732426   15556105.1     15600740
      NMVSDIQEATK                           5370097    5897534.2      6212194
      VLGIDAGEEK                                 NA           NA           NA

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
      00125amol_R1 00125amol_R1
      00125amol_R2 00125amol_R2
      ...                   ...
      00500amol_R2 00500amol_R2
      00500amol_R3 00500amol_R3
      00050amol_R1 00050amol_R1
      00050amol_R2 00050amol_R2
      00050amol_R3 00050amol_R3

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
      colnames(27): 12500amol_R1 12500amol_R2 ... 00050amol_R2 00050amol_R3
      colData names(1): sample

---

    Code
      SummarizedExperiment::assays(D_norm_lts)$intensities_norm
    Output
                                       12500amol_R1 12500amol_R2 12500amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                       8121984    8345342.1    5604632.4
      AEWALR                                5251358    5112729.3    5435498.5
      DEGLHTDFACLLFAHLK                     4696671    3914991.2     528532.7
      DIHDWNNR                              2001959    1977337.8           NA
      ELETLREENR                                 NA           NA           NA
      ESEFLFNAIHTIPEIGEK                   34550852   39049894.8   36017247.5
      GMMPGLTFSNELICR                       9209961    7932258.6    8895769.7
      IVTEAVEIEQR                          19220262   21125229.5   20888913.8
      LLVAFGNK                             12057575   13809961.4   12898642.7
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                           17022103   18075229.4   18864101.5
      TNFFEK                                7256143    7792341.4    5430376.5
      TVLFPIK                              14371798   14706893.7   14000410.9
      VENPFDFMENISLAGK                     17382407   20041167.4   21214758.2
      WIQDADALFGER                         18275598   20806893.8   20061770.3
      YFLDALPVALLGMNADLMNQYVEFVADR         26878492   32918706.4   30967204.4
      AANLGGVAVSGLEMAQNSQK                  9619517    3840999.7   10837976.7
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           4520052           NA    5206426.6
      FLGFEQIFK                            41897223   43773308.1   43903989.4
      GANIASFVMVADAMLDQGDVF                58346069   59909486.0   61495227.3
      GCIISETGITSEQIHDIASAK                 8698571    8985024.8    9121899.2
      GGLCVDLK                                   NA   11210578.4           NA
      ICYAFMR                               9684311   10332466.6   10558992.8
      NSWEGVLTGK                           12594499   13073272.4   15663416.0
      SLEEIVDEYSTFSESK                      8981768    9165594.3   11150634.7
      VLPIVSVPER                           41394614   40156541.5   41737505.7
      VTISGSGNVAQYAALK                      4999549     148208.9    5547963.9
      VTWENDNGEQEVAQGYR                     3700536           NA    3444142.5
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 29946629   34485649.3   30003838.4
      ALVAQGVK                             11437892   14408991.7   13567768.0
      FIAEGSNMGSTPEAIAVFETAR               24224150   44430413.2   46827870.7
      GANIASFIK                            19751130   17263258.4   18525179.8
      GCIISETGITSEQVADISSAK                 3440754    4043401.0    4085805.3
      HIGQDTDVPAGDIGVGGR                   28859661   34984662.0   36274435.7
      IMINCFNECIDYAK                       14645306   15128473.4   13649501.5
      ITWTSER                              20902690   19216290.8   18099075.6
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        11204754   94984486.5   99612481.9
      SLEQIVNEYSTFSENK                     55865319   74562606.8   61039699.0
      STATGPSEAVWYGPPK                     50854366   43252786.2   42673627.2
      VDIALPCATQNEVSGEEAK                  50027785   49164581.1   53839518.8
      VIELGGTVVSLSDSK                      14160864   16110366.9   17500786.0
      VQYIAGARPWTHVQK                            NA    7285155.3    4518339.2
      VTWENDKGEQEVAQGYR                    35172554   35511636.6   40672790.1
      AAGLTAAYAR                           65641977   65506171.6   56502942.5
      APEAEQVLSAAATFPIAQPATDVEAR           29984981   31998114.0   30957396.4
      AVQDNGESAFR                          10971616   14472443.8   13583024.9
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        28330811   30911901.0   30548728.7
      GFTLAEVK                             29421816   32107810.7   31115414.6
      IAPRPLDLLRPVVR                       17893091   18521544.7   20674227.0
      IIVFPR                               52116944   52123167.2   53231421.3
      NQEIFDANVQR                          87870829   94618830.7   90919302.4
      TIGIAVDHR                            25354516   28857775.4   25774400.0
      VHFDQAGK                              9630215   10492817.5   12941144.2
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           168252983  189140870.1  191409057.7
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR            113177910  119988894.5  128583208.4
      CCSDVFNQVVK                          51762696   53594394.3   51291611.9
      DIVGAVLK                                   NA   56130325.2           NA
      EALDFFAR                            193756869  203046547.3  201227980.5
      EKDIVGAVLK                           39747509   41097567.6   44300669.6
      GVIFYESHGK                           90579671  102004003.5  109533844.2
      IGDYAGIK                            144787792  139476200.8  148940309.2
      LPLVGGHEGAGVVVGMGENVK               123119077  131904973.6  154748839.6
      SANLMAGHWVAISGAAGGLGSLAVQYAK         28994901   32243318.5   25549905.2
      SIGGEVFIDFTK                         73284060   71273855.0   61501765.9
      SIPETQK                              30700543   28829813.5   28671036.8
      SISIVGSYVGNR                        265827803  300321760.2  280999917.6
      VLGIDGGEGK                           22178389   22527627.3           NA
      VLGIDGGEGKEELFR                     297205556  325476732.3  341025030.9
      VVGLSTLPEIYEK                       543272137  542375183.9  563296674.0
      YSGVCHTDLHAWHGDWPLPVK               101187955   99891803.2  111092230.4
      ANGTVVLVGLPAGAK                       5263772    4864728.6    4709268.7
      CSSDVFNHVVK                                NA           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        40763829   39281118.3   40314252.3
      DIPVPEPKPNEILINVK                    19052725   23058903.7   18569860.8
      EALDFFSR                              7674480    6355206.4    6701060.7
      GVIFYENK                              5675044    5731332.9    5464050.8
      IQQGTDLAEVAPILCAGVTVYK                8421026    8727990.2    6343503.7
      IVGLSELPK                            17274417   18019305.6   20267738.9
      NMVSDIQEATK                           5916760    7379903.1    6027576.3
      VLGIDAGEEK                                 NA           NA           NA
                                       00125amol_R1 00125amol_R2 00125amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                      13218450   14668357.4   13788629.4
      AEWALR                                6207258    6424893.7    5476021.2
      DEGLHTDFACLLFAHLK                     3077795    5325002.9    4292601.1
      DIHDWNNR                              2126158           NA    2625754.2
      ELETLREENR                            5318855    5325422.5    6854984.7
      ESEFLFNAIHTIPEIGEK                   34560611   34920614.3   41209976.4
      GMMPGLTFSNELICR                       6759773    8298117.3    8229321.2
      IVTEAVEIEQR                          18004130   20813479.2   19350150.5
      LLVAFGNK                              6662884    8153982.8    7955419.5
      LLVAFGNKK                             1706084    2228629.4    1877478.6
      NKPDPAIVEK                           19822461   20341423.1   17749792.4
      TNFFEK                                5933181    5660477.5    5130802.5
      TVLFPIK                              11876033   11427954.3   10808201.2
      VENPFDFMENISLAGK                     19185382   18986097.5   15865775.5
      WIQDADALFGER                         17080744   18002122.7   17703521.8
      YFLDALPVALLGMNADLMNQYVEFVADR         25939182   24538526.5   21947337.6
      AANLGGVAVSGLEMAQNSQK                 10161986   10379989.7   10110119.2
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                           25282194   28075800.5   24609907.1
      FHPSVNLSILK                           5020603    5295315.8    5266194.2
      FLGFEQIFK                            50978679   49844931.0   48901955.5
      GANIASFVMVADAMLDQGDVF                82452488   91526438.1   82439064.6
      GCIISETGITSEQIHDIASAK                 8556773    8338714.1    8986951.4
      GGLCVDLK                             11219045           NA           NA
      ICYAFMR                              12005914   14651573.2   12307971.2
      NSWEGVLTGK                           13496224   19407800.9           NA
      SLEEIVDEYSTFSESK                      9566620   10805889.3    9368784.2
      VLPIVSVPER                           34406081   42320356.4   33621603.2
      VTISGSGNVAQYAALK                      6961989     125220.8           NA
      VTWENDNGEQEVAQGYR                          NA    3270090.1    2436849.5
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA    1911678.6
      AANLGGVAVSGLEMAQNSQR                  9922133   38491456.7   35513667.1
      ALVAQGVK                             14857602   15426794.3   14331805.6
      FIAEGSNMGSTPEAIAVFETAR               42101280   46106246.5   44989074.8
      GANIASFIK                            19170213   20602627.4   20223256.0
      GCIISETGITSEQVADISSAK                 3007167    3399643.3    2827735.3
      HIGQDTDVPAGDIGVGGR                   35657487   38815601.9   34521867.5
      IMINCFNECIDYAK                       12909391   14141752.6   12862212.1
      ITWTSER                              17884678   18682932.5     550911.5
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR         9749591  122619201.7  102992249.9
      SLEQIVNEYSTFSENK                     47646336   56047748.6   46676944.7
      STATGPSEAVWYGPPK                     39445835   47229740.1   48589126.2
      VDIALPCATQNEVSGEEAK                  48477763   57223692.8   51690260.2
      VIELGGTVVSLSDSK                      13512341   14419741.2   13641770.6
      VQYIAGARPWTHVQK                       3799819           NA    3900206.6
      VTWENDKGEQEVAQGYR                    44446718   45035203.6   46200156.6
      AAGLTAAYAR                           60414470   65885398.3   63966043.1
      APEAEQVLSAAATFPIAQPATDVEAR           20326816   26831670.3   21973490.6
      AVQDNGESAFR                          14902160   14849836.8   15155622.9
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        41710690   47044064.7   42789211.1
      GFTLAEVK                             32933783   36923181.3   33192091.6
      IAPRPLDLLRPVVR                       23824116   30963734.9   28780293.5
      IIVFPR                               42209356   48123499.7    3359746.2
      NQEIFDANVQR                          81425766   86907631.2   82135288.3
      TIGIAVDHR                            20481345   21278192.2   18073686.4
      VHFDQAGK                                   NA           NA   14186958.6
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           164758979  188958822.9  171593397.0
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             82960635   68544647.8   83175370.2
      CCSDVFNQVVK                          57960007   64652807.3   56210693.5
      DIVGAVLK                                   NA   62858994.0   52919447.9
      EALDFFAR                            198110845  218173852.1  203892265.5
      EKDIVGAVLK                           32617139   36758486.2   35871761.1
      GVIFYESHGK                          117783851  126143887.5  111190187.6
      IGDYAGIK                            135550025  144323293.5  141507468.7
      LPLVGGHEGAGVVVGMGENVK               113185882   90758560.1  152089346.5
      SANLMAGHWVAISGAAGGLGSLAVQYAK         23233490   31786161.6   33891179.6
      SIGGEVFIDFTK                         81313898   67515565.4   78098684.1
      SIPETQK                              27784058   51081718.1   49200702.4
      SISIVGSYVGNR                        252148833  297888396.9  279041703.8
      VLGIDGGEGK                           10133545   13010915.9   11214577.5
      VLGIDGGEGKEELFR                     403929214  412503624.8  331719740.0
      VVGLSTLPEIYEK                       528036458  579758355.9  545308707.2
      YSGVCHTDLHAWHGDWPLPVK                91661699  112831904.7  119126595.9
      ANGTVVLVGLPAGAK                       4787860    6578993.8    5169931.3
      CSSDVFNHVVK                           4185669    5506482.3    4351847.5
      DIPVPKPKPNELLINVK                          NA    1378403.9           NA
      VVGLSSLPEIYEK                        36980944   39750273.0   36808237.1
      DIPVPEPKPNEILINVK                    18493316   23719246.9   21442384.9
      EALDFFSR                              6922266    8048766.8    6661553.6
      GVIFYENK                                   NA    6644976.8    5866705.7
      IQQGTDLAEVAPILCAGVTVYK                8018858    8335672.0    8397806.4
      IVGLSELPK                            19105747   21319103.7   21847755.3
      NMVSDIQEATK                           3736205    4790215.8           NA
      VLGIDAGEEK                                 NA           NA           NA
                                       25000amol_R1 25000amol_R2 25000amol_R3
      AAADALSDLEIK                               NA           NA     34122959
      AAADALSDLEIKDSK                       7953482   10317170.1      5866018
      AEWALR                                5103988    5812357.5      5903302
      DEGLHTDFACLLFAHLK                     5402317    3521704.6      6225033
      DIHDWNNR                              2726314     685679.2      2162948
      ELETLREENR                                 NA    1143065.6           NA
      ESEFLFNAIHTIPEIGEK                   34760064   31695718.3     38055768
      GMMPGLTFSNELICR                       8475171    9276786.4      8048825
      IVTEAVEIEQR                          20341205   23435011.1     25115463
      LLVAFGNK                             11600325   16088078.1     13094112
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                           16367173   17908090.3     21699566
      TNFFEK                                6079805    7051793.3      8074524
      TVLFPIK                              13007933   15813003.8     14602014
      VENPFDFMENISLAGK                     19062969   18021134.5     16332839
      WIQDADALFGER                         20132220   19921219.6     19007928
      YFLDALPVALLGMNADLMNQYVEFVADR         32694540   24873499.0     24743599
      AANLGGVAVSGLEMAQNSQK                  9715922   11152755.4      9292527
      DAVWFGPPK                           176697360           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           6259045    7700196.2      7225855
      FLGFEQIFK                            39319033   42030787.1     36837960
      GANIASFVMVADAMLDQGDVF                56560847   47203502.7     48099033
      GCIISETGITSEQIHDIASAK                 8177173    7326019.8      9431441
      GGLCVDLK                             10832831   12534721.1     12941278
      ICYAFMR                               9153985    7214577.0      9868041
      NSWEGVLTGK                           14841251   16261412.6     14496879
      SLEEIVDEYSTFSESK                     10676147    7759732.8      9623700
      VLPIVSVPER                           41815795   55830661.6     44091267
      VTISGSGNVAQYAALK                      1996304    5677740.7      1051248
      VTWENDNGEQEVAQGYR                     2666936           NA      2069787
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 32999725   34811028.9     30145370
      ALVAQGVK                             13275523   18854835.7     18802526
      FIAEGSNMGSTPEAIAVFETAR               40851589   47196908.5     40021199
      GANIASFIK                            18871676   19132736.1           NA
      GCIISETGITSEQVADISSAK                 3204326    3367493.4      3617406
      HIGQDTDVPAGDIGVGGR                   31459428   37529742.7     38539581
      IMINCFNECIDYAK                       12223963   12838056.5     12288081
      ITWTSER                              18755573   23798636.7     22538695
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR       104055726   76991599.6     90674611
      SLEQIVNEYSTFSENK                     58334455   55956894.3     60935176
      STATGPSEAVWYGPPK                     44131216   55586674.5     44255783
      VDIALPCATQNEVSGEEAK                  52634362   59336916.9     53322660
      VIELGGTVVSLSDSK                      14813608   17334390.8     16032038
      VQYIAGARPWTHVQK                       3820555    6991503.1      3103707
      VTWENDKGEQEVAQGYR                    35795037   33634426.9     33962337
      AAGLTAAYAR                           60934051   74514988.9     73746010
      APEAEQVLSAAATFPIAQPATDVEAR           33394474   32855363.7     35411830
      AVQDNGESAFR                          12509244   14929374.8     15766281
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        20285918   26121695.7     25566179
      GFTLAEVK                             31554522   31463035.6     30019793
      IAPRPLDLLRPVVR                       18219289   15502132.2     20529458
      IIVFPR                               53966779   65742756.6     55992882
      NQEIFDANVQR                          85252606   98762976.5     94541224
      TIGIAVDHR                            26568737   28101853.8     30892020
      VHFDQAGK                             12628664   11908267.7     13234292
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           188075432  207040509.1    189037669
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             87632160   97387605.0    119989613
      CCSDVFNQVVK                          48847200   56222548.3     51710598
      DIVGAVLK                             51053153   63331146.3     54670913
      EALDFFAR                            201034712  233690686.5    190079279
      EKDIVGAVLK                           40216894   35103059.8     34120039
      GVIFYESHGK                           99587428  101956476.0    103489338
      IGDYAGIK                            151353784  189650538.2    180315399
      LPLVGGHEGAGVVVGMGENVK               130112505  105620993.2    124048973
      SANLMAGHWVAISGAAGGLGSLAVQYAK         29685820   48419670.2     35185012
      SIGGEVFIDFTK                         96653898   63397088.8     53802579
      SIPETQK                              25384489   23521678.4     15469373
      SISIVGSYVGNR                        288034389  312557878.4    292556211
      VLGIDGGEGK                           21648191   31837965.6     32111775
      VLGIDGGEGKEELFR                     320532107  293792536.0    291592478
      VVGLSTLPEIYEK                       529307004  593821345.4    532885886
      YSGVCHTDLHAWHGDWPLPVK               118104162   96530352.9    114664746
      ANGTVVLVGLPAGAK                       4668769    5489145.2      5111678
      CSSDVFNHVVK                           4542825           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        41929686   46440454.2     41589455
      DIPVPEPKPNEILINVK                    22910947   26036912.5     24573242
      EALDFFSR                              6752427    7511412.3      5593837
      GVIFYENK                              4905176    4624639.5      6610040
      IQQGTDLAEVAPILCAGVTVYK                8612946    7781776.5      6805122
      IVGLSELPK                            20112317   19565130.3     18163932
      NMVSDIQEATK                           6114082    7185091.3      6761608
      VLGIDAGEEK                                 NA           NA           NA
                                       02500amol_R1 02500amol_R2 02500amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                     6252233.7    8474606.1      8946931
      AEWALR                              6070218.4    5284894.3      5016018
      DEGLHTDFACLLFAHLK                   5367006.8    5702891.6      3942675
      DIHDWNNR                             683491.6           NA           NA
      ELETLREENR                                 NA           NA      4779228
      ESEFLFNAIHTIPEIGEK                 37374713.7   29614298.8     37005204
      GMMPGLTFSNELICR                     9086000.5   10881900.9      9664109
      IVTEAVEIEQR                        19100979.5   20655069.9     19176343
      LLVAFGNK                           11694000.5   11933706.8      9894192
      LLVAFGNKK                            828201.1           NA           NA
      NKPDPAIVEK                         15649699.0   16584792.5     18762606
      TNFFEK                              7227937.2    6496855.5      6357103
      TVLFPIK                            11057902.9   13110334.5     12727826
      VENPFDFMENISLAGK                   18140992.6   17923076.6     16833210
      WIQDADALFGER                       20834053.3   24613342.4     19891355
      YFLDALPVALLGMNADLMNQYVEFVADR       26872201.4   26712566.3     25659941
      AANLGGVAVSGLEMAQNSQK               11685505.0   11851723.9      8703641
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                         5905300.3    5376923.0      8064877
      FLGFEQIFK                          43136759.3   47888400.8     40930030
      GANIASFVMVADAMLDQGDVF              54503684.6   63443211.1     59678192
      GCIISETGITSEQIHDIASAK               6609892.5    9126774.3      8887192
      GGLCVDLK                           10586086.3   10719667.2     11707412
      ICYAFMR                             9182317.8   10009302.8     10120768
      NSWEGVLTGK                         13793441.0   15213022.5     12261469
      SLEEIVDEYSTFSESK                    9104584.3   12397891.6      8619346
      VLPIVSVPER                         44483289.7   44767277.5     44648503
      VTISGSGNVAQYAALK                    7442341.4     197867.4           NA
      VTWENDNGEQEVAQGYR                   3658888.4    3252757.8      3307109
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR               34578645.4   38515793.3     34187460
      ALVAQGVK                                   NA   14162255.8     13942211
      FIAEGSNMGSTPEAIAVFETAR             39169379.5   42332039.2     37671723
      GANIASFIK                          19435488.3   19993433.3     19152613
      GCIISETGITSEQVADISSAK               4039060.2    3606900.8      3725180
      HIGQDTDVPAGDIGVGGR                 36930825.9   32213502.9     33445623
      IMINCFNECIDYAK                     14073791.2   15574440.0     14164040
      ITWTSER                                    NA   19247504.5     17386235
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR      95615338.6  107817894.4     87880396
      SLEQIVNEYSTFSENK                   56757105.3   59326746.2     61366155
      STATGPSEAVWYGPPK                   52933086.4   56182529.1     46753299
      VDIALPCATQNEVSGEEAK                57721340.0   63880838.0     56617983
      VIELGGTVVSLSDSK                    14962628.7   16125226.5     14205310
      VQYIAGARPWTHVQK                     8894534.1           NA      6569956
      VTWENDKGEQEVAQGYR                  36795960.4   36379619.4     35243985
      AAGLTAAYAR                         59032826.6   62290831.8     58583491
      APEAEQVLSAAATFPIAQPATDVEAR         28631823.5   31747008.7     30101679
      AVQDNGESAFR                        14357327.1   14791561.1     14136182
      DGKAPEAEQVLSAAATFPIAQPATDVEAR      30450913.9   25616766.8     21977579
      GFTLAEVK                           32414427.0   35864628.3     31290270
      IAPRPLDLLRPVVR                      9479977.4   22366318.4     10578664
      IIVFPR                             63949574.0   63787308.3     53628966
      NQEIFDANVQR                        92314852.5   97582734.5     90831237
      TIGIAVDHR                          24075071.3   26064786.0     24250553
      VHFDQAGK                            9906661.9           NA     11276136
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                         165523415.7  163746392.0    160263714
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR           76760727.4  103595198.1     93652077
      CCSDVFNQVVK                        56174104.4   58016174.6     55658444
      DIVGAVLK                           51081076.2           NA     51319881
      EALDFFAR                          208722829.0  221388447.4    205124760
      EKDIVGAVLK                         31887708.5   34634884.9     37846091
      GVIFYESHGK                         91446191.7   96966130.8     95441153
      IGDYAGIK                          146971455.3  144024310.5    148584159
      LPLVGGHEGAGVVVGMGENVK             108625073.4  129151921.2    159799421
      SANLMAGHWVAISGAAGGLGSLAVQYAK       40473432.6   24060246.6     29850961
      SIGGEVFIDFTK                       56072158.8   79989898.9     76279255
      SIPETQK                            43250386.1   24239223.3     24836594
      SISIVGSYVGNR                      280955476.4  293764329.6    266896399
      VLGIDGGEGK                         20904140.9   22034922.8     22030199
      VLGIDGGEGKEELFR                   297999492.7  342850142.5    345011159
      VVGLSTLPEIYEK                     565075511.5  579053244.8    523846604
      YSGVCHTDLHAWHGDWPLPVK             157527064.2  169566022.6    148480983
      ANGTVVLVGLPAGAK                     5357555.6    4455366.0      5318221
      CSSDVFNHVVK                                NA           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                      39342474.5   43207293.1     36082808
      DIPVPEPKPNEILINVK                  23606759.0   26726422.6     23038231
      EALDFFSR                            7725665.0    7233546.8      9221999
      GVIFYENK                            5358829.9    5848266.9      5918191
      IQQGTDLAEVAPILCAGVTVYK              7199265.1    9527797.6      7275166
      IVGLSELPK                          20685382.8   22046469.7     20501127
      NMVSDIQEATK                         6761536.5    7681681.5      6387540
      VLGIDAGEEK                                 NA           NA           NA
                                       00250amol_R1 00250amol_R2 00250amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                    11025078.7     14410611   13471469.8
      AEWALR                              5324150.3      6410555    5476389.6
      DEGLHTDFACLLFAHLK                          NA           NA           NA
      DIHDWNNR                                   NA      1913884    2340844.1
      ELETLREENR                          4275483.1      4197643    4417215.1
      ESEFLFNAIHTIPEIGEK                 28762650.8     25587626   39910244.1
      GMMPGLTFSNELICR                     7870974.3           NA    8248658.3
      IVTEAVEIEQR                        18815451.9     17958059   20834207.0
      LLVAFGNK                            7233202.2      6734732    7848744.1
      LLVAFGNKK                                  NA      1824942           NA
      NKPDPAIVEK                         17118137.8     20709119   20293724.7
      TNFFEK                              6934786.2      5984043    5880576.4
      TVLFPIK                             9706461.8     10308044   12522955.3
      VENPFDFMENISLAGK                   18249325.2     18813740   19610110.6
      WIQDADALFGER                       15067927.4     15650684   18498169.6
      YFLDALPVALLGMNADLMNQYVEFVADR       24846674.2     24307683   26203781.9
      AANLGGVAVSGLEMAQNSQK               11148752.3     10898945    6426506.4
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA     21794822   23459712.3
      FHPSVNLSILK                         5010381.7      6876732    5806233.4
      FLGFEQIFK                          48504588.9     47222967   47792100.7
      GANIASFVMVADAMLDQGDVF              77710759.0     87501320   83419075.4
      GCIISETGITSEQIHDIASAK               9002481.6      9345224    7555858.2
      GGLCVDLK                           11874801.9     12075634   11516760.8
      ICYAFMR                            10753209.9     11375345   13125390.2
      NSWEGVLTGK                         12737318.9     14667213   15331113.7
      SLEEIVDEYSTFSESK                    9057068.5      8684704    8836246.0
      VLPIVSVPER                         34545459.7     34806901   42749378.8
      VTISGSGNVAQYAALK                    3441432.3           NA     172954.4
      VTWENDNGEQEVAQGYR                   3046956.0      2816489    3471370.9
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR               38011520.1     36195210   42635087.1
      ALVAQGVK                           14383457.7     13675564   13102959.1
      FIAEGSNMGSTPEAIAVFETAR             44085388.5     28660718   46149290.7
      GANIASFIK                          18353808.1     19058074   20685734.6
      GCIISETGITSEQVADISSAK               4310239.6      3834515    3801214.7
      HIGQDTDVPAGDIGVGGR                 34413256.8     38250072   33986300.9
      IMINCFNECIDYAK                     15028479.7     14008840   14206355.0
      ITWTSER                            18113923.8     17841514   17935256.1
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     108054513.3    102487497  105462424.8
      SLEQIVNEYSTFSENK                   55879162.8     55793069   64365469.4
      STATGPSEAVWYGPPK                   48125038.8     47051217   41072388.0
      VDIALPCATQNEVSGEEAK                50379950.4     44890035   54465883.1
      VIELGGTVVSLSDSK                    14201145.7     13891273   14223445.3
      VQYIAGARPWTHVQK                    11026144.8     10574869    4935693.6
      VTWENDKGEQEVAQGYR                  37723659.1     43867716   42238804.5
      AAGLTAAYAR                         64746349.7     59801582   63270618.8
      APEAEQVLSAAATFPIAQPATDVEAR         20840074.7     19263560   16115133.6
      AVQDNGESAFR                          286453.7     14537379   12610543.4
      DGKAPEAEQVLSAAATFPIAQPATDVEAR      28312734.6     28258946   31391771.9
      GFTLAEVK                           31734016.7     35858867   33736354.5
      IAPRPLDLLRPVVR                     13439913.1     25521176   26020060.7
      IIVFPR                             44904193.3     50141688   47092464.5
      NQEIFDANVQR                        93233933.3     85388186   83934990.4
      TIGIAVDHR                          21808140.8     24116509   22360589.1
      VHFDQAGK                           11227647.6     11446908   11164272.3
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                         168953120.2    175695772  170006274.8
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR           79985927.5     79958649   70722012.2
      CCSDVFNQVVK                        53333191.7     59148320   59895274.3
      DIVGAVLK                           55180833.2     57496252           NA
      EALDFFAR                          196470505.4    190222928  218350607.6
      EKDIVGAVLK                         35114784.9     38251094   40084352.1
      GVIFYESHGK                        102289829.0     94612572  107487631.5
      IGDYAGIK                          138215955.7    134608765           NA
      LPLVGGHEGAGVVVGMGENVK             122276982.1    126440435   91813642.6
      SANLMAGHWVAISGAAGGLGSLAVQYAK       27732747.9     31829907   31803008.5
      SIGGEVFIDFTK                       65806104.9     63070959   57491943.6
      SIPETQK                            28859670.7     24754436   27119183.9
      SISIVGSYVGNR                      252368858.9    261304779  277664809.9
      VLGIDGGEGK                         10425048.3     11393747   11463353.5
      VLGIDGGEGKEELFR                   325549537.7    338949920  370241109.5
      VVGLSTLPEIYEK                     548545906.8    509329633  565284749.7
      YSGVCHTDLHAWHGDWPLPVK             101948660.4    116810188  117528213.2
      ANGTVVLVGLPAGAK                     4902807.0      4872782    5797901.9
      CSSDVFNHVVK                                NA      4033151    4433557.7
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                      38145855.3     35375310   36621419.5
      DIPVPEPKPNEILINVK                    786553.7     20506700   22737645.0
      EALDFFSR                            7646229.5           NA    7494760.2
      GVIFYENK                                   NA      5865250    6537059.6
      IQQGTDLAEVAPILCAGVTVYK              9972573.3      6734323    9931203.4
      IVGLSELPK                          18471084.8     19161328   21528502.6
      NMVSDIQEATK                         6141995.5      5737971    6963570.7
      VLGIDAGEEK                                 NA           NA           NA
                                       50000amol_R1 50000amol_R2 50000amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                       4182266      2898907     11060373
      AEWALR                                4843557      5146672      5378537
      DEGLHTDFACLLFAHLK                          NA      5038796      2646551
      DIHDWNNR                              2677100      1280767           NA
      ELETLREENR                            5166613      5172267      6069318
      ESEFLFNAIHTIPEIGEK                   37423416     37506351     32842839
      GMMPGLTFSNELICR                       8756814      8170045      6471505
      IVTEAVEIEQR                          24636300     25961341     24317635
      LLVAFGNK                             15448247     16848756     14988059
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                           21278519     24414297     24050472
      TNFFEK                                5304714      6257094      5990803
      TVLFPIK                              15503523     15754829     15151432
      VENPFDFMENISLAGK                     13173169     14064645     14353786
      WIQDADALFGER                         14853375     13862733     17537643
      YFLDALPVALLGMNADLMNQYVEFVADR         25043410     25603967     25878331
      AANLGGVAVSGLEMAQNSQK                  9760112      9707894      9241643
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                          12511615      6319468      8262749
      FLGFEQIFK                            38441191     42847067     40826988
      GANIASFVMVADAMLDQGDVF                42397353     44464260     47019797
      GCIISETGITSEQIHDIASAK                 7440550      8285978     10000368
      GGLCVDLK                             11690376     14168918     12250115
      ICYAFMR                               9152080      9491763      9326982
      NSWEGVLTGK                           15749193     15690369     14762219
      SLEEIVDEYSTFSESK                      5501689      7577201      5957840
      VLPIVSVPER                           44330248     44823530     39827528
      VTISGSGNVAQYAALK                      5881775      5189993           NA
      VTWENDNGEQEVAQGYR                          NA      2676046           NA
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 30794717     24359316     26015757
      ALVAQGVK                             15664963     15734922     15292702
      FIAEGSNMGSTPEAIAVFETAR               19342992     31295422     35058950
      GANIASFIK                            19736941     20366576     20313068
      GCIISETGITSEQVADISSAK                 3648197      4104597      5096286
      HIGQDTDVPAGDIGVGGR                   42950110     47851263     38607033
      IMINCFNECIDYAK                       15283297     11973480     11633140
      ITWTSER                              16416888     18468792           NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        75361861     10850167     10358828
      SLEQIVNEYSTFSENK                     54045614     50445786     58159934
      STATGPSEAVWYGPPK                     40829453           NA     42879726
      VDIALPCATQNEVSGEEAK                  54031576     53972138     47775158
      VIELGGTVVSLSDSK                      14226917     13615319     14200984
      VQYIAGARPWTHVQK                       4157348      4663884           NA
      VTWENDKGEQEVAQGYR                    34182329     41423255     42091690
      AAGLTAAYAR                           64778754     69586671     71885214
      APEAEQVLSAAATFPIAQPATDVEAR           28449447     32383662     28862296
      AVQDNGESAFR                          13821562     13857993     15972143
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        14316411     12588696     15554099
      GFTLAEVK                             31471186     31157022     35087780
      IAPRPLDLLRPVVR                       19921193     11775360     12716209
      IIVFPR                               57945623     59082504     54050615
      NQEIFDANVQR                          82651237     92756317     94014608
      TIGIAVDHR                                  NA     26172733     29298599
      VHFDQAGK                             19212260     24308127     23593988
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           223884206    227127374    223869474
      ANGTTVLVGMPAGAK                      21744415     19990242     21361540
      ATDGGAHGVINVSVSEAAIEASTR            101356364    131498787    163594336
      CCSDVFNQVVK                          49964862     56788479     53729635
      DIVGAVLK                             56066248     58901446     59626450
      EALDFFAR                            196860521    204689541    201593044
      EKDIVGAVLK                           41574359     46192362     42447267
      GVIFYESHGK                          112701047    129337475    144854457
      IGDYAGIK                            149323404    160581708    162594876
      LPLVGGHEGAGVVVGMGENVK               139751053    153851306    126566258
      SANLMAGHWVAISGAAGGLGSLAVQYAK         19660608     19829092     38058291
      SIGGEVFIDFTK                         55413469     65234660     49511721
      SIPETQK                              13286352     10549668      9864864
      SISIVGSYVGNR                        305297442    319021062    308179708
      VLGIDGGEGK                           23250897     26339571     29363948
      VLGIDGGEGKEELFR                     330689176    347430943    329466288
      VVGLSTLPEIYEK                       519188148    483261834    541198087
      YSGVCHTDLHAWHGDWPLPVK               105550300     99590993    113678987
      ANGTVVLVGLPAGAK                       5849136      6363832      4207055
      CSSDVFNHVVK                           4266496      3940224      4227621
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        34886875     38009709     38768484
      DIPVPEPKPNEILINVK                    21088125     20905008     19831597
      EALDFFSR                              7179614      7688868           NA
      GVIFYENK                                   NA      7260398      4923302
      IQQGTDLAEVAPILCAGVTVYK                6868841      7342964      6604222
      IVGLSELPK                            21618947     20197842     21290424
      NMVSDIQEATK                           5338669      4447469      6591536
      VLGIDAGEEK                                 NA           NA           NA
                                       05000amol_R1 05000amol_R2 05000amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                            NA      2289399           NA
      AEWALR                                5848626      5048690      5464988
      DEGLHTDFACLLFAHLK                     5425231           NA           NA
      DIHDWNNR                                   NA           NA           NA
      ELETLREENR                                 NA      4564929           NA
      ESEFLFNAIHTIPEIGEK                   29330027     33466564     36643775
      GMMPGLTFSNELICR                      12035649     10330899      9409919
      IVTEAVEIEQR                          16531288     17190632     19350869
      LLVAFGNK                             10605354           NA     11917068
      LLVAFGNKK                                  NA      1432927           NA
      NKPDPAIVEK                           15989934     16959153     18852490
      TNFFEK                                6480772      7664505      6519640
      TVLFPIK                              11768967     12452960     11462726
      VENPFDFMENISLAGK                     19181129     17095797     19595763
      WIQDADALFGER                         19614612     21591792     21423870
      YFLDALPVALLGMNADLMNQYVEFVADR         23447042     23486453     25254087
      AANLGGVAVSGLEMAQNSQK                  9160777      6026408     10319354
      DAVWFGPPK                           163784706           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           5541093      4950898      5294100
      FLGFEQIFK                            41582418     44610188     41301045
      GANIASFVMVADAMLDQGDVF                54327221     51644296     57489783
      GCIISETGITSEQIHDIASAK                 7514540      7904958      7223494
      GGLCVDLK                             11174675           NA      9590044
      ICYAFMR                               9761061     10162031      9620441
      NSWEGVLTGK                           10805116     12618156     13498134
      SLEEIVDEYSTFSESK                      9045714      8510678      9307343
      VLPIVSVPER                           44117395     46242778     45688717
      VTISGSGNVAQYAALK                           NA      6648853      3188876
      VTWENDNGEQEVAQGYR                     3792678      3381531      3802506
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 36434556     31527035      9562333
      ALVAQGVK                             11861856     13512460     13228536
      FIAEGSNMGSTPEAIAVFETAR               40970148     39914326     40433177
      GANIASFIK                            19303982     16911225     20097364
      GCIISETGITSEQVADISSAK                 4431716      4241369      2870837
      HIGQDTDVPAGDIGVGGR                   33182434     32107262     29604166
      IMINCFNECIDYAK                       15083015     12071580     14334852
      ITWTSER                              17838729     21789620     19494798
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        81868381     94918636     84837279
      SLEQIVNEYSTFSENK                     63129725     59669581     58912957
      STATGPSEAVWYGPPK                     44609808     49293815     53337337
      VDIALPCATQNEVSGEEAK                  53132645     50111639     59057960
      VIELGGTVVSLSDSK                      13851484     14800381     14139367
      VQYIAGARPWTHVQK                      13949367      4420026     10147714
      VTWENDKGEQEVAQGYR                    34821479     36467633     37547088
      AAGLTAAYAR                           56606503     61041120     60781880
      APEAEQVLSAAATFPIAQPATDVEAR           28562942     28103797     29849060
      AVQDNGESAFR                          14342898     14622948     13279019
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        32357417     28567774     29317383
      GFTLAEVK                             28590909     30994735     30715853
      IAPRPLDLLRPVVR                       19422840     17262013     10246423
      IIVFPR                               54981441     55402006     59185777
      NQEIFDANVQR                          97198102     90684711    100089407
      TIGIAVDHR                            28342205     28688103     26944711
      VHFDQAGK                              9468810      9590370     10076179
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           164923348    170733624    173852781
      ANGTTVLVGMPAGAK                      22404285     18677400     21804099
      ATDGGAHGVINVSVSEAAIEASTR             83873990     91492338    103718769
      CCSDVFNQVVK                          49478003     51403639     53575786
      DIVGAVLK                                   NA     55391808           NA
      EALDFFAR                            198553248    205506464    206891899
      EKDIVGAVLK                           31750142     34842182     37336566
      GVIFYESHGK                           89052815     88828800    100630750
      IGDYAGIK                            138994265    148727821    140931814
      LPLVGGHEGAGVVVGMGENVK               126479188    131790083    138450658
      SANLMAGHWVAISGAAGGLGSLAVQYAK         25289845     15024722     24650447
      SIGGEVFIDFTK                         66306937     66413065     63614264
      SIPETQK                              35019244     34220146     19522724
      SISIVGSYVGNR                        272435152    280283366    291756030
      VLGIDGGEGK                           19490760     21272616     22087660
      VLGIDGGEGKEELFR                     300182067    311160832    323989577
      VVGLSTLPEIYEK                       541903842    546932923    583544283
      YSGVCHTDLHAWHGDWPLPVK               123652559    139213727    159664006
      ANGTVVLVGLPAGAK                       4169628      4985161      4512374
      CSSDVFNHVVK                                NA      3670115      1422422
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        40292955     41813065     41225859
      DIPVPEPKPNEILINVK                    19898274     21666233     22068326
      EALDFFSR                              6800391      6923262      7012112
      GVIFYENK                                   NA           NA           NA
      IQQGTDLAEVAPILCAGVTVYK                7534716      9245904      8595541
      IVGLSELPK                            18578847     20549627     20127439
      NMVSDIQEATK                           5482663      6017639      7575689
      VLGIDAGEEK                                 NA           NA      4698300
                                       00500amol_R1 00500amol_R2 00500amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                      12894341      3729110    3625605.2
      AEWALR                                5679358      5376759    5481807.8
      DEGLHTDFACLLFAHLK                          NA           NA           NA
      DIHDWNNR                              1961258           NA           NA
      ELETLREENR                            4314729      3482088    6463946.5
      ESEFLFNAIHTIPEIGEK                   32026500     34877920   36028414.5
      GMMPGLTFSNELICR                      10443837      8333672    8701601.4
      IVTEAVEIEQR                          18909041     17532258   19286424.1
      LLVAFGNK                              7453387      8867548   10041364.8
      LLVAFGNKK                             2311791           NA           NA
      NKPDPAIVEK                           18670732     20216142   22674392.8
      TNFFEK                                6932239      5483882    6762109.0
      TVLFPIK                               9635864      9836395   13510389.5
      VENPFDFMENISLAGK                     17678567     15007899   18915183.5
      WIQDADALFGER                         16966488     15350152   16660083.5
      YFLDALPVALLGMNADLMNQYVEFVADR         27221386     25581009   26007047.7
      AANLGGVAVSGLEMAQNSQK                 10571062      9698141   11470162.2
      DAVWFGPPK                                  NA           NA  184662909.9
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           4794576      5233187    5208323.8
      FLGFEQIFK                            42422865     47894177   48749519.0
      GANIASFVMVADAMLDQGDVF                69597729     69278233   77506299.6
      GCIISETGITSEQIHDIASAK                 6679783      9156046    9151344.9
      GGLCVDLK                             10469472     10089701           NA
      ICYAFMR                              11454041     12897338   12683874.1
      NSWEGVLTGK                           13548505     12440033   16383514.6
      SLEEIVDEYSTFSESK                     10004246      8296159    9440678.5
      VLPIVSVPER                           40553418     44389465   45073280.4
      VTISGSGNVAQYAALK                           NA           NA     122200.9
      VTWENDNGEQEVAQGYR                     3170749      2568543    2929715.4
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 33857020      9476256   40110996.6
      ALVAQGVK                             12793700     14422975   15395312.8
      FIAEGSNMGSTPEAIAVFETAR               39223253     42583451   45643437.8
      GANIASFIK                            20490808     21303811   24358271.7
      GCIISETGITSEQVADISSAK                 4026289      3443028    3755805.3
      HIGQDTDVPAGDIGVGGR                   34549162     41913447   36864503.5
      IMINCFNECIDYAK                       13110813     13658223   15075131.1
      ITWTSER                              17827629     19443655   21965951.0
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        93360306     97281119  110414804.4
      SLEQIVNEYSTFSENK                     54290862     50973489   56036045.1
      STATGPSEAVWYGPPK                     43400788     48537110   48859082.8
      VDIALPCATQNEVSGEEAK                  57339512     55048622   59263391.1
      VIELGGTVVSLSDSK                      14854934     15471971   14852812.3
      VQYIAGARPWTHVQK                       6215507      9024172           NA
      VTWENDKGEQEVAQGYR                    39500490     40538634   42916043.2
      AAGLTAAYAR                           64237192     66351679   72809947.3
      APEAEQVLSAAATFPIAQPATDVEAR           19275525     18156822   15476156.0
      AVQDNGESAFR                          12625649     14387203   16447338.2
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        23816695     28624549   28894001.3
      GFTLAEVK                             30520689     34461222   37877171.1
      IAPRPLDLLRPVVR                       29924441     24103713   31181013.1
      IIVFPR                               48257172     50683444   53842641.2
      NQEIFDANVQR                          93888194     89319381   95945997.4
      TIGIAVDHR                            22655530     22627383   26700597.4
      VHFDQAGK                             11405619     11326260   11724392.8
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           173823422    164387881  161750242.4
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             82479726     93512466  107329997.8
      CCSDVFNQVVK                          48121402     54906500   56882771.4
      DIVGAVLK                             58306042     48053702   54513852.6
      EALDFFAR                            204832123    210708080  222914512.7
      EKDIVGAVLK                           36578114     36877298   37065547.8
      GVIFYESHGK                           96472569     98944044  104995118.2
      IGDYAGIK                            125240592    123597875  144858266.6
      LPLVGGHEGAGVVVGMGENVK               141703877    119334213  161112006.5
      SANLMAGHWVAISGAAGGLGSLAVQYAK         38790309     29971325   39908888.6
      SIGGEVFIDFTK                         55389363     68537652   85288522.2
      SIPETQK                              24252488     20573864   21314950.4
      SISIVGSYVGNR                        273533580    280376908  296513745.1
      VLGIDGGEGK                           10191286      8591425   10657475.2
      VLGIDGGEGKEELFR                     324091245    355382564  387653826.4
      VVGLSTLPEIYEK                       523378640    584421493  590634105.6
      YSGVCHTDLHAWHGDWPLPVK               123227780     89276841  126764279.8
      ANGTVVLVGLPAGAK                       4914300      5676279    5348735.7
      CSSDVFNHVVK                           4394861           NA    1140953.0
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        35167247     33132816   38115445.8
      DIPVPEPKPNEILINVK                          NA     19746269   23960438.1
      EALDFFSR                              6652724      6674131   10971274.5
      GVIFYENK                              5217361      5019810    6584041.2
      IQQGTDLAEVAPILCAGVTVYK                8740922      8416335    8624694.0
      IVGLSELPK                            19158744     22145908   23547712.2
      NMVSDIQEATK                           6185695      6663882    5222045.8
      VLGIDAGEEK                                 NA           NA           NA
                                       00050amol_R1 00050amol_R2 00050amol_R3
      AAADALSDLEIK                         16458044           NA           NA
      AAADALSDLEIKDSK                       1745328           NA           NA
      AEWALR                                4528760    5588844.9      5232325
      DEGLHTDFACLLFAHLK                    10922386    7824804.0      5850027
      DIHDWNNR                                   NA    2057646.2      2023635
      ELETLREENR                            6697556    7744526.6      7260462
      ESEFLFNAIHTIPEIGEK                   30216421   33442814.7     27362483
      GMMPGLTFSNELICR                       4847661           NA           NA
      IVTEAVEIEQR                          24127591   23349098.6     21236910
      LLVAFGNK                              6404343    6552082.2      7262483
      LLVAFGNKK                             1861951    2565306.9      3262991
      NKPDPAIVEK                           17866190   16701727.0     20444913
      TNFFEK                                4906051    6720051.4      5637327
      TVLFPIK                               8375589    9882267.1     10769137
      VENPFDFMENISLAGK                     11595333   14064564.0     14572927
      WIQDADALFGER                         13592110   13549214.2     14637242
      YFLDALPVALLGMNADLMNQYVEFVADR         26047144   29115890.2     27366158
      AANLGGVAVSGLEMAQNSQK                       NA   10810445.9      8253766
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                           25474273   23486403.1     26445530
      FHPSVNLSILK                          11130416    7870297.5      7894335
      FLGFEQIFK                            51796427   56498996.8     53646306
      GANIASFVMVADAMLDQGDVF               104251616   90894705.1     85397848
      GCIISETGITSEQIHDIASAK                 8493788   10467184.5      7709291
      GGLCVDLK                              9923997   10393040.0     11064068
      ICYAFMR                               8727822   10010418.0      9775006
      NSWEGVLTGK                           10244711   13756086.5     14259620
      SLEEIVDEYSTFSESK                      8798742    7825261.6      5921693
      VLPIVSVPER                           28712928   31126944.4     33444872
      VTISGSGNVAQYAALK                      3522570    3055850.3           NA
      VTWENDNGEQEVAQGYR                     2827717    2893007.1      3181770
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK      7534011    5658321.0      4798197
      AANLGGVAVSGLEMAQNSQR                 25750858   27722706.6     29069594
      ALVAQGVK                                   NA     692362.8           NA
      FIAEGSNMGSTPEAIAVFETAR               34872678   39989496.2     39298483
      GANIASFIK                            14052298   17509993.2     17376155
      GCIISETGITSEQVADISSAK                 6411041    5107546.7           NA
      HIGQDTDVPAGDIGVGGR                   38005742   36458938.2     33477948
      IMINCFNECIDYAK                       11061861   14055410.4     12993527
      ITWTSER                                    NA   18510485.8     16984751
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        79863200           NA     97391732
      SLEQIVNEYSTFSENK                     42403537    9553651.5     44893913
      STATGPSEAVWYGPPK                     36527465   38242066.8     37987369
      VDIALPCATQNEVSGEEAK                   7098566   44605675.6     45836592
      VIELGGTVVSLSDSK                      16465136   15127301.4     13455678
      VQYIAGARPWTHVQK                            NA           NA           NA
      VTWENDKGEQEVAQGYR                    51260591   50184817.7     48708729
      AAGLTAAYAR                           43710031   56878872.8     54219631
      APEAEQVLSAAATFPIAQPATDVEAR            4192286    2753505.7      2893729
      AVQDNGESAFR                          14286333   16903107.0     14230218
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        26448233   23005837.2     29570335
      GFTLAEVK                             27382794   34436899.7     33692945
      IAPRPLDLLRPVVR                       16731478   10017740.9      9707016
      IIVFPR                                1766762    2017919.4     42535379
      NQEIFDANVQR                          64031611   78713044.3     78044772
      TIGIAVDHR                            27867411   29863742.4     25228134
      VHFDQAGK                             11088653   14557945.1     14589465
      VHFDQAGKK                             2235618    1919792.4      1662734
      ANELLINVK                           145723745  163740267.8    167973174
      ANGTTVLVGMPAGAK                      17870918   18842762.8     19157689
      ATDGGAHGVINVSVSEAAIEASTR             98144662   92488353.4     83228585
      CCSDVFNQVVK                          56587431   60969633.4     63463576
      DIVGAVLK                             47312740   55936963.5     53940319
      EALDFFAR                            178937695  192601686.9    195086664
      EKDIVGAVLK                           49192106   32457883.3     31166273
      GVIFYESHGK                          180600363  207741803.3    157195768
      IGDYAGIK                                   NA  118750139.4    118156384
      LPLVGGHEGAGVVVGMGENVK               100886882  100735780.8    125341321
      SANLMAGHWVAISGAAGGLGSLAVQYAK         48631842   42175384.8     40876964
      SIGGEVFIDFTK                         70603481   40619266.5     59844463
      SIPETQK                              26484481   29417960.2     54712102
      SISIVGSYVGNR                        236035746  243038229.6    249626034
      VLGIDGGEGK                            8885421    8714354.4      9044017
      VLGIDGGEGKEELFR                     357702031  421579929.3    448332731
      VVGLSTLPEIYEK                       443995263  482350928.7    511931369
      YSGVCHTDLHAWHGDWPLPVK                66942461  113642409.7     97501987
      ANGTVVLVGLPAGAK                       4956719    5377670.5      5159005
      CSSDVFNHVVK                           5165616    4065954.2      5037909
      DIPVPKPKPNELLINVK                     6754685    6036640.9      4124999
      VVGLSSLPEIYEK                        26969885   29260517.7     32925756
      DIPVPEPKPNEILINVK                     3354018   25649407.7     19338690
      EALDFFSR                              6168418    6381000.7      5635306
      GVIFYENK                              5752751    6954384.6      6225169
      IQQGTDLAEVAPILCAGVTVYK                6102147    7862700.0      7791614
      IVGLSELPK                            13120101   16046326.6     15663531
      NMVSDIQEATK                           5998369    7102490.4      6914262
      VLGIDAGEEK                                 NA           NA           NA

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
      00125amol_R1 00125amol_R1
      00125amol_R2 00125amol_R2
      ...                   ...
      00500amol_R2 00500amol_R2
      00500amol_R3 00500amol_R3
      00050amol_R1 00050amol_R1
      00050amol_R2 00050amol_R2
      00050amol_R3 00050amol_R3

