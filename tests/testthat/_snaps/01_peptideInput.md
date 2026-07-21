# read MaxQuant Output table

    Code
      D1
    Output
      class: SummarizedExperiment 
      dim: 87 27 
      metadata(0):
      assays(1): intensities
      rownames(87): AAADALSDLEIK AAADALSDLEIKDSK ... NMVSDIQEATK VLGIDAGEEK
      rowData names(1): Sequence
      colnames(27): 12500amol_R1 12500amol_R2 ... 00050amol_R2 00050amol_R3
      colData names(1): sample

---

    Code
      SummarizedExperiment::assays(D1)$intensities
    Output
                                       12500amol_R1 12500amol_R2 12500amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                       8047500      7759800      5142900
      AEWALR                                5203200      4754000      4987700
      DEGLHTDFACLLFAHLK                     4653600      3640300       484990
      DIHDWNNR                              1983600      1838600           NA
      ELETLREENR                                 NA           NA           NA
      ESEFLFNAIHTIPEIGEK                   34234000     36310000     33050000
      GMMPGLTFSNELICR                       9125500      7375700      8162900
      IVTEAVEIEQR                          19044000     19643000     19168000
      LLVAFGNK                             11947000     12841000     11836000
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                           16866000     16807000     17310000
      TNFFEK                                7189600      7245600      4983000
      TVLFPIK                              14240000     13675000     12847000
      VENPFDFMENISLAGK                     17223000     18635000     19467000
      WIQDADALFGER                         18108000     19347000     18409000
      YFLDALPVALLGMNADLMNQYVEFVADR         26632000     30609000     28416000
      AANLGGVAVSGLEMAQNSQK                  9531300      3571500      9945100
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           4478600           NA      4777500
      FLGFEQIFK                            41513000     40702000     40287000
      GANIASFVMVADAMLDQGDVF                57811000     55706000     56429000
      GCIISETGITSEQIHDIASAK                 8618800      8354600      8370400
      GGLCVDLK                                   NA     10424000           NA
      ICYAFMR                               9595500      9607500      9689100
      NSWEGVLTGK                           12479000     12156000     14373000
      SLEEIVDEYSTFSESK                      8899400      8522500     10232000
      VLPIVSVPER                           41015000     37339000     38299000
      VTISGSGNVAQYAALK                      4953700       137810      5090900
      VTWENDNGEQEVAQGYR                     3666600           NA      3160400
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 29672000     32066000     27532000
      ALVAQGVK                             11333000     13398000     12450000
      FIAEGSNMGSTPEAIAVFETAR               24002000     41313000     42970000
      GANIASFIK                            19570000     16052000     16999000
      GCIISETGITSEQVADISSAK                 3409200      3759700      3749200
      HIGQDTDVPAGDIGVGGR                   28595000     32530000     33286000
      IMINCFNECIDYAK                       14511000     14067000     12525000
      ITWTSER                              20711000     17868000     16608000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        11102000     88320000     91406000
      SLEQIVNEYSTFSENK                     55353000     69331000     56011000
      STATGPSEAVWYGPPK                     50388000     40218000     39158000
      VDIALPCATQNEVSGEEAK                  49569000     45715000     49404000
      VIELGGTVVSLSDSK                      14031000     14980000     16059000
      VQYIAGARPWTHVQK                            NA      6774000      4146100
      VTWENDKGEQEVAQGYR                    34850000     33020000     37322000
      AAGLTAAYAR                           65040000     60910000     51848000
      APEAEQVLSAAATFPIAQPATDVEAR           29710000     29753000     28407000
      AVQDNGESAFR                          10871000     13457000     12464000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        28071000     28743000     28032000
      GFTLAEVK                             29152000     29855000     28552000
      IAPRPLDLLRPVVR                       17729000     17222000     18971000
      IIVFPR                               51639000     48466000     48846000
      NQEIFDANVQR                          87065000     87980000     83429000
      TIGIAVDHR                            25122000     26833000     23651000
      VHFDQAGK                              9541900      9756600     11875000
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           166710000    175870000    175640000
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR            112140000    111570000    117990000
      CCSDVFNQVVK                          51288000     49834000     47066000
      DIVGAVLK                                   NA     52192000           NA
      EALDFFAR                            191980000    188800000    184650000
      EKDIVGAVLK                           39383000     38214000     40651000
      GVIFYESHGK                           89749000     94847000    100510000
      IGDYAGIK                            143460000    129690000    136670000
      LPLVGGHEGAGVVVGMGENVK               121990000    122650000    142000000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         28729000     29981000     23445000
      SIGGEVFIDFTK                         72612000     66273000     56435000
      SIPETQK                              30419000     26807000     26309000
      SISIVGSYVGNR                        263390000    279250000    257850000
      VLGIDGGEGK                           21975000     20947000           NA
      VLGIDGGEGKEELFR                     294480000    302640000    312930000
      VVGLSTLPEIYEK                       538290000    504320000    516890000
      YSGVCHTDLHAWHGDWPLPVK               100260000     92883000    101940000
      ANGTVVLVGLPAGAK                       5215500      4523400      4321300
      CSSDVFNHVVK                                NA           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        40390000     36525000     36993000
      DIPVPEPKPNEILINVK                    18878000     21441000     17040000
      EALDFFSR                              7604100      5909300      6149000
      GVIFYENK                              5623000      5329200      5013900
      IQQGTDLAEVAPILCAGVTVYK                8343800      8115600      5820900
      IVGLSELPK                            17116000     16755000     18598000
      NMVSDIQEATK                           5862500      6862100      5531000
      VLGIDAGEEK                                 NA           NA           NA
                                       00125amol_R1 00125amol_R2 00125amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                      13943000     13983000     13708000
      AEWALR                                6547500      6124700      5444000
      DEGLHTDFACLLFAHLK                     3246500      5076200      4267500
      DIHDWNNR                              2242700           NA      2610400
      ELETLREENR                            5610400      5076600      6814900
      ESEFLFNAIHTIPEIGEK                   36455000     33289000     40969000
      GMMPGLTFSNELICR                       7130300      7910400      8181200
      IVTEAVEIEQR                          18991000     19841000     19237000
      LLVAFGNK                              7028100      7773000      7908900
      LLVAFGNKK                             1799600      2124500      1866500
      NKPDPAIVEK                           20909000     19391000     17646000
      TNFFEK                                6258400      5396000      5100800
      TVLFPIK                              12527000     10894000     10745000
      VENPFDFMENISLAGK                     20237000     18099000     15773000
      WIQDADALFGER                         18017000     17161000     17600000
      YFLDALPVALLGMNADLMNQYVEFVADR         27361000     23392000     21819000
      AANLGGVAVSGLEMAQNSQK                 10719000      9895000     10051000
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                           26668000     26764000     24466000
      FHPSVNLSILK                           5295800      5047900      5235400
      FLGFEQIFK                            53773000     47516000     48616000
      GANIASFVMVADAMLDQGDVF                86972000     87250000     81957000
      GCIISETGITSEQIHDIASAK                 9025800      7949100      8934400
      GGLCVDLK                             11834000           NA           NA
      ICYAFMR                              12664000     13967000     12236000
      NSWEGVLTGK                           14236000     18501000           NA
      SLEEIVDEYSTFSESK                     10091000     10301000      9314000
      VLPIVSVPER                           36292000     40343000     33425000
      VTISGSGNVAQYAALK                      7343600       119370           NA
      VTWENDNGEQEVAQGYR                          NA      3117300      2422600
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA      1900500
      AANLGGVAVSGLEMAQNSQR                 10466000     36693000     35306000
      ALVAQGVK                             15672000     14706000     14248000
      FIAEGSNMGSTPEAIAVFETAR               44409000     43952000     44726000
      GANIASFIK                            20221000     19640000     20105000
      GCIISETGITSEQVADISSAK                 3172000      3240800      2811200
      HIGQDTDVPAGDIGVGGR                   37612000     37002000     34320000
      IMINCFNECIDYAK                       13617000     13481000     12787000
      ITWTSER                              18865000     17810000       547690
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        10284000    116890000    102390000
      SLEQIVNEYSTFSENK                     50258000     53429000     46404000
      STATGPSEAVWYGPPK                     41608000     45023000     48305000
      VDIALPCATQNEVSGEEAK                  51135000     54550000     51388000
      VIELGGTVVSLSDSK                      14253000     13746000     13562000
      VQYIAGARPWTHVQK                       4008100           NA      3877400
      VTWENDKGEQEVAQGYR                    46883000     42931000     45930000
      AAGLTAAYAR                           63726000     62807000     63592000
      APEAEQVLSAAATFPIAQPATDVEAR           21441000     25578000     21845000
      AVQDNGESAFR                          15719000     14156000     15067000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        43997000     44846000     42539000
      GFTLAEVK                             34739000     35198000     32998000
      IAPRPLDLLRPVVR                       25130000     29517000     28612000
      IIVFPR                               44523000     45875000      3340100
      NQEIFDANVQR                          85889000     82847000     81655000
      TIGIAVDHR                            21604000     20284000     17968000
      VHFDQAGK                                   NA           NA     14104000
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           173790000    180130000    170590000
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             87508000     65342000     82689000
      CCSDVFNQVVK                          61137000     61632000     55882000
      DIVGAVLK                                   NA     59922000     52610000
      EALDFFAR                            208970000    207980000    202700000
      EKDIVGAVLK                           34405000     35041000     35662000
      GVIFYESHGK                          124240000    120250000    110540000
      IGDYAGIK                            142980000    137580000    140680000
      LPLVGGHEGAGVVVGMGENVK               119390000     86518000    151200000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         24507000     30301000     33693000
      SIGGEVFIDFTK                         85771000     64361000     77642000
      SIPETQK                              29307000     48695000     48913000
      SISIVGSYVGNR                        265970000    283970000    277410000
      VLGIDGGEGK                           10689000     12403000     11149000
      VLGIDGGEGKEELFR                     426070000    393230000    329780000
      VVGLSTLPEIYEK                       556980000    552670000    542120000
      YSGVCHTDLHAWHGDWPLPVK                96686000    107560000    118430000
      ANGTVVLVGLPAGAK                       5050300      6271600      5139700
      CSSDVFNHVVK                           4415100      5249200      4326400
      DIPVPKPKPNELLINVK                          NA      1314000           NA
      VVGLSSLPEIYEK                        39008000     37893000     36593000
      DIPVPEPKPNEILINVK                    19507000     22611000     21317000
      EALDFFSR                              7301700      7672700      6622600
      GVIFYENK                                   NA      6334500      5832400
      IQQGTDLAEVAPILCAGVTVYK                8458400      7946200      8348700
      IVGLSELPK                            20153000     20323000     21720000
      NMVSDIQEATK                           3941000      4566400           NA
      VLGIDAGEEK                                 NA           NA           NA
                                       25000amol_R1 25000amol_R2 25000amol_R3
      AAADALSDLEIK                               NA           NA     35053000
      AAADALSDLEIKDSK                       7192900     10952000      6025900
      AEWALR                                4615900      6170000      6064200
      DEGLHTDFACLLFAHLK                     4885700      3738400      6394700
      DIHDWNNR                              2465600       727870      2221900
      ELETLREENR                                 NA      1213400           NA
      ESEFLFNAIHTIPEIGEK                   31436000     33646000     39093000
      GMMPGLTFSNELICR                       7664700      9847600      8268200
      IVTEAVEIEQR                          18396000     24877000     25800000
      LLVAFGNK                             10491000     17078000     13451000
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                           14802000     19010000     22291000
      TNFFEK                                5498400      7485700      8294600
      TVLFPIK                              11764000     16786000     15000000
      VENPFDFMENISLAGK                     17240000     19130000     16778000
      WIQDADALFGER                         18207000     21147000     19526000
      YFLDALPVALLGMNADLMNQYVEFVADR         29568000     26404000     25418000
      AANLGGVAVSGLEMAQNSQK                  8786800     11839000      9545800
      DAVWFGPPK                           159800000           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           5660500      8174000      7422800
      FLGFEQIFK                            35559000     44617000     37842000
      GANIASFVMVADAMLDQGDVF                51152000     50108000     49410000
      GCIISETGITSEQIHDIASAK                 7395200      7776800      9688500
      GGLCVDLK                              9796900     13306000     13294000
      ICYAFMR                               8278600      7658500     10137000
      NSWEGVLTGK                           13422000     17262000     14892000
      SLEEIVDEYSTFSESK                      9655200      8237200      9886000
      VLPIVSVPER                           37817000     59266000     45293000
      VTISGSGNVAQYAALK                      1805400      6027100      1079900
      VTWENDNGEQEVAQGYR                     2411900           NA      2126200
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 29844000     36953000     30967000
      ALVAQGVK                             12006000     20015000     19315000
      FIAEGSNMGSTPEAIAVFETAR               36945000     50101000     41112000
      GANIASFIK                            17067000     20310000           NA
      GCIISETGITSEQVADISSAK                 2897900      3574700      3716000
      HIGQDTDVPAGDIGVGGR                   28451000     39839000     39590000
      IMINCFNECIDYAK                       11055000     13628000     12623000
      ITWTSER                              16962000     25263000     23153000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        94105000     81729000     93146000
      SLEQIVNEYSTFSENK                     52756000     59400000     62596000
      STATGPSEAVWYGPPK                     39911000     59007000     45462000
      VDIALPCATQNEVSGEEAK                  47601000     62988000     54776000
      VIELGGTVVSLSDSK                      13397000     18401000     16469000
      VQYIAGARPWTHVQK                       3455200      7421700      3188300
      VTWENDKGEQEVAQGYR                    32372000     35704000     34888000
      AAGLTAAYAR                           55107000     79100000     75756000
      APEAEQVLSAAATFPIAQPATDVEAR           30201000     34877000     36377000
      AVQDNGESAFR                          11313000     15848000     16196000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        18346000     27729000     26263000
      GFTLAEVK                             28537000     33399000     30838000
      IAPRPLDLLRPVVR                       16477000     16456000     21089000
      IIVFPR                               48806000     69788000     57519000
      NQEIFDANVQR                          77100000    104840000     97118000
      TIGIAVDHR                            24028000     29831000     31734000
      VHFDQAGK                             11421000     12641000     13595000
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           170090000    219780000    194190000
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             79252000    103380000    123260000
      CCSDVFNQVVK                          44176000     59682000     53120000
      DIVGAVLK                             46171000     67228000     56161000
      EALDFFAR                            181810000    248070000    195260000
      EKDIVGAVLK                           36371000     37263000     35050000
      GVIFYESHGK                           90064000    108230000    106310000
      IGDYAGIK                            136880000    201320000    185230000
      LPLVGGHEGAGVVVGMGENVK               117670000    112120000    127430000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         26847000     51399000     36144000
      SIGGEVFIDFTK                         87411000     67298000     55269000
      SIPETQK                              22957000     24969000     15891000
      SISIVGSYVGNR                        260490000    331790000    300530000
      VLGIDGGEGK                           19578000     33797000     32987000
      VLGIDGGEGKEELFR                     289880000    311870000    299540000
      VVGLSTLPEIYEK                       478690000    630360000    547410000
      YSGVCHTDLHAWHGDWPLPVK               106810000    102470000    117790000
      ANGTVVLVGLPAGAK                       4222300      5826900      5251000
      CSSDVFNHVVK                           4108400           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        37920000     49298000     42723000
      DIPVPEPKPNEILINVK                    20720000     27639000     25243000
      EALDFFSR                              6106700      7973600      5746300
      GVIFYENK                              4436100      4909200      6790200
      IQQGTDLAEVAPILCAGVTVYK                7789300      8260600      6990600
      IVGLSELPK                            18189000     20769000     18659000
      NMVSDIQEATK                           5529400      7627200      6945900
      VLGIDAGEEK                                 NA           NA           NA
                                       02500amol_R1 02500amol_R2 02500amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                       5887600      7339300      8671500
      AEWALR                                5716200      4576900      4861600
      DEGLHTDFACLLFAHLK                     5054000      4938900      3821300
      DIHDWNNR                               643630           NA           NA
      ELETLREENR                                 NA           NA      4632100
      ESEFLFNAIHTIPEIGEK                   35195000     25647000     35866000
      GMMPGLTFSNELICR                       8556100      9424100      9366600
      IVTEAVEIEQR                          17987000     17888000     18586000
      LLVAFGNK                             11012000     10335000      9589600
      LLVAFGNKK                              779900           NA           NA
      NKPDPAIVEK                           14737000     14363000     18185000
      TNFFEK                                6806400      5626500      6161400
      TVLFPIK                              10413000     11354000     12336000
      VENPFDFMENISLAGK                     17083000     15522000     16315000
      WIQDADALFGER                         19619000     21316000     19279000
      YFLDALPVALLGMNADLMNQYVEFVADR         25305000     23134000     24870000
      AANLGGVAVSGLEMAQNSQK                 11004000     10264000      8435700
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           5560900      4656600      7816600
      FLGFEQIFK                            40621000     41473000     39670000
      GANIASFVMVADAMLDQGDVF                51325000     54944000     57841000
      GCIISETGITSEQIHDIASAK                 6224400      7904100      8613600
      GGLCVDLK                              9968700      9283600     11347000
      ICYAFMR                               8646800      8668400      9809200
      NSWEGVLTGK                           12989000     13175000     11884000
      SLEEIVDEYSTFSESK                      8573600     10737000      8354000
      VLPIVSVPER                           41889000     38770000     43274000
      VTISGSGNVAQYAALK                      7008300       171360           NA
      VTWENDNGEQEVAQGYR                     3445500      2817000      3205300
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 32562000     33356000     33135000
      ALVAQGVK                                   NA     12265000     13513000
      FIAEGSNMGSTPEAIAVFETAR               36885000     36661000     36512000
      GANIASFIK                            18302000     17315000     18563000
      GCIISETGITSEQVADISSAK                 3803500      3123700      3610500
      HIGQDTDVPAGDIGVGGR                   34777000     27898000     32416000
      IMINCFNECIDYAK                       13253000     13488000     13728000
      ITWTSER                                    NA     16669000     16851000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        90039000     93374000     85175000
      SLEQIVNEYSTFSENK                     53447000     51379000     59477000
      STATGPSEAVWYGPPK                     49846000     48656000     45314000
      VDIALPCATQNEVSGEEAK                  54355000     55323000     54875000
      VIELGGTVVSLSDSK                      14090000     13965000     13768000
      VQYIAGARPWTHVQK                       8375800           NA      6367700
      VTWENDKGEQEVAQGYR                    34650000     31506000     34159000
      AAGLTAAYAR                           55590000     53946000     56780000
      APEAEQVLSAAATFPIAQPATDVEAR           26962000     27494000     29175000
      AVQDNGESAFR                          13520000     12810000     13701000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        28675000     22185000     21301000
      GFTLAEVK                             30524000     31060000     30327000
      IAPRPLDLLRPVVR                        8927100     19370000     10253000
      IIVFPR                               60220000     55242000     51978000
      NQEIFDANVQR                          86931000     84510000     88035000
      TIGIAVDHR                            22671000     22573000     23504000
      VHFDQAGK                              9328900           NA     10929000
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           155870000    141810000    155330000
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             72284000     89717000     90769000
      CCSDVFNQVVK                          52898000     50244000     53945000
      DIVGAVLK                             48102000           NA     49740000
      EALDFFAR                            196550000    191730000    198810000
      EKDIVGAVLK                           30028000     29995000     36681000
      GVIFYESHGK                           86113000     83976000     92503000
      IGDYAGIK                            138400000    124730000    144010000
      LPLVGGHEGAGVVVGMGENVK               102290000    111850000    154880000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         38113000     20837000     28932000
      SIGGEVFIDFTK                         52802000     69274000     73931000
      SIPETQK                              40728000     20992000     24072000
      SISIVGSYVGNR                        264570000    254410000    258680000
      VLGIDGGEGK                           19685000     19083000     21352000
      VLGIDGGEGKEELFR                     280620000    296920000    334390000
      VVGLSTLPEIYEK                       532120000    501480000    507720000
      YSGVCHTDLHAWHGDWPLPVK               148340000    146850000    143910000
      ANGTVVLVGLPAGAK                       5045100      3858500      5154500
      CSSDVFNHVVK                                NA           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        37048000     37419000     34972000
      DIPVPEPKPNEILINVK                    22230000     23146000     22329000
      EALDFFSR                              7275100      6264500      8938100
      GVIFYENK                              5046300      5064800      5736000
      IQQGTDLAEVAPILCAGVTVYK                6779400      8251400      7051200
      IVGLSELPK                            19479000     19093000     19870000
      NMVSDIQEATK                           6367200      6652600      6190900
      VLGIDAGEEK                                 NA           NA           NA
                                       00250amol_R1 00250amol_R2 00250amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                      10341000     14096000     12612000
      AEWALR                                4993800      6270600      5127000
      DEGLHTDFACLLFAHLK                          NA           NA           NA
      DIHDWNNR                                   NA      1872100      2191500
      ELETLREENR                            4010200      4106000      4135400
      ESEFLFNAIHTIPEIGEK                   26978000     25029000     37364000
      GMMPGLTFSNELICR                       7382600           NA      7722400
      IVTEAVEIEQR                          17648000     17566000     19505000
      LLVAFGNK                              6784400      6587700      7348000
      LLVAFGNKK                                  NA      1785100           NA
      NKPDPAIVEK                           16056000     20257000     18999000
      TNFFEK                                6504500      5853400      5505400
      TVLFPIK                               9104200     10083000     11724000
      VENPFDFMENISLAGK                     17117000     18403000     18359000
      WIQDADALFGER                         14133000     15309000     17318000
      YFLDALPVALLGMNADLMNQYVEFVADR         23305000     23777000     24532000
      AANLGGVAVSGLEMAQNSQK                 10457000     10661000      6016500
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA     21319000     21963000
      FHPSVNLSILK                           4699500      6726600      5435800
      FLGFEQIFK                            45495000     46192000     44743000
      GANIASFVMVADAMLDQGDVF                72889000     85591000     78097000
      GCIISETGITSEQIHDIASAK                 8443900      9141200      7073800
      GGLCVDLK                             11138000     11812000     10782000
      ICYAFMR                              10086000     11127000     12288000
      NSWEGVLTGK                           11947000     14347000     14353000
      SLEEIVDEYSTFSESK                      8495100      8495100      8272500
      VLPIVSVPER                           32402000     34047000     40022000
      VTISGSGNVAQYAALK                      3227900           NA       161920
      VTWENDNGEQEVAQGYR                     2857900      2755000      3249900
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 35653000     35405000     39915000
      ALVAQGVK                             13491000     13377000     12267000
      FIAEGSNMGSTPEAIAVFETAR               41350000     28035000     43205000
      GANIASFIK                            17215000     18642000     19366000
      GCIISETGITSEQVADISSAK                 4042800      3750800      3558700
      HIGQDTDVPAGDIGVGGR                   32278000     37415000     31818000
      IMINCFNECIDYAK                       14096000     13703000     13300000
      ITWTSER                              16990000     17452000     16791000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR       101350000    100250000     98734000
      SLEQIVNEYSTFSENK                     52412000     54575000     60259000
      STATGPSEAVWYGPPK                     45139000     46024000     38452000
      VDIALPCATQNEVSGEEAK                  47254000     43910000     50991000
      VIELGGTVVSLSDSK                      13320000     13588000     13316000
      VQYIAGARPWTHVQK                      10342000     10344000      4620800
      VTWENDKGEQEVAQGYR                    35383000     42910000     39544000
      AAGLTAAYAR                           60729000     58496000     59234000
      APEAEQVLSAAATFPIAQPATDVEAR           19547000     18843000     15087000
      AVQDNGESAFR                            268680     14220000     11806000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        26556000     27642000     29389000
      GFTLAEVK                             29765000     35076000     31584000
      IAPRPLDLLRPVVR                       12606000     24964000     24360000
      IIVFPR                               42118000     49047000     44088000
      NQEIFDANVQR                          87449000     83524000     78580000
      TIGIAVDHR                            20455000     23590000     20934000
      VHFDQAGK                             10531000     11197000     10452000
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           158470000    171860000    159160000
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             75023000     78213000     66210000
      CCSDVFNQVVK                          50024000     57857000     56074000
      DIVGAVLK                             51757000     56241000           NA
      EALDFFAR                            184280000    186070000    204420000
      EKDIVGAVLK                           32936000     37416000     37527000
      GVIFYESHGK                           95943000     92547000    100630000
      IGDYAGIK                            129640000    131670000           NA
      LPLVGGHEGAGVVVGMGENVK               114690000    123680000     85956000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         26012000     31135000     29774000
      SIGGEVFIDFTK                         61723000     61694000     53824000
      SIPETQK                              27069000     24214000     25389000
      SISIVGSYVGNR                        236710000    255600000    259950000
      VLGIDGGEGK                            9778200     11145000     10732000
      VLGIDGGEGKEELFR                     305350000    331550000    346620000
      VVGLSTLPEIYEK                       514510000    498210000    529220000
      YSGVCHTDLHAWHGDWPLPVK                95623000    114260000    110030000
      ANGTVVLVGLPAGAK                       4598600      4766400      5428000
      CSSDVFNHVVK                                NA      3945100      4150700
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        35779000     34603000     34285000
      DIPVPEPKPNEILINVK                      737750     20059000     21287000
      EALDFFSR                              7171800           NA      7016600
      GVIFYENK                                   NA      5737200      6120000
      IQQGTDLAEVAPILCAGVTVYK                9353800      6587300      9297600
      IVGLSELPK                            17325000     18743000     20155000
      NMVSDIQEATK                           5760900      5612700      6519300
      VLGIDAGEEK                                 NA           NA           NA
                                       50000amol_R1 50000amol_R2 50000amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                       4766700      3058100     11509000
      AEWALR                                5520400      5429300      5596700
      DEGLHTDFACLLFAHLK                          NA      5315500      2753900
      DIHDWNNR                              3051200      1351100           NA
      ELETLREENR                            5888600      5456300      6315500
      ESEFLFNAIHTIPEIGEK                   42653000     39566000     34175000
      GMMPGLTFSNELICR                       9980500      8618700      6734000
      IVTEAVEIEQR                          28079000     27387000     25304000
      LLVAFGNK                             17607000     17774000     15596000
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                           24252000     25755000     25026000
      TNFFEK                                6046000      6600700      6233800
      TVLFPIK                              17670000     16620000     15766000
      VENPFDFMENISLAGK                     15014000     14837000     14936000
      WIQDADALFGER                         16929000     14624000     18249000
      YFLDALPVALLGMNADLMNQYVEFVADR         28543000     27010000     26928000
      AANLGGVAVSGLEMAQNSQK                 11124000     10241000      9616500
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                          14260000      6666500      8597900
      FLGFEQIFK                            43813000     45200000     42483000
      GANIASFVMVADAMLDQGDVF                48322000     46906000     48927000
      GCIISETGITSEQIHDIASAK                 8480300      8741000     10406000
      GGLCVDLK                             13324000     14947000     12747000
      ICYAFMR                              10431000     10013000      9705300
      NSWEGVLTGK                           17950000     16552000     15361000
      SLEEIVDEYSTFSESK                      6270500      7993300      6199500
      VLPIVSVPER                           50525000     47285000     41443000
      VTISGSGNVAQYAALK                      6703700      5475000           NA
      VTWENDNGEQEVAQGYR                          NA      2823000           NA
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 35098000     25697000     27071000
      ALVAQGVK                             17854000     16599000     15913000
      FIAEGSNMGSTPEAIAVFETAR               22046000     33014000     36481000
      GANIASFIK                            22495000     21485000     21137000
      GCIISETGITSEQVADISSAK                 4158000      4330000      5303000
      HIGQDTDVPAGDIGVGGR                   48952000     50479000     40173000
      IMINCFNECIDYAK                       17419000     12631000     12105000
      ITWTSER                              18711000     19483000           NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        85893000     11446000     10779000
      SLEQIVNEYSTFSENK                     61598000     53216000     60519000
      STATGPSEAVWYGPPK                     46535000           NA     44619000
      VDIALPCATQNEVSGEEAK                  61582000     56936000     49713000
      VIELGGTVVSLSDSK                      16215000     14363000     14777000
      VQYIAGARPWTHVQK                       4738300      4920000           NA
      VTWENDKGEQEVAQGYR                    38959000     43698000     43799000
      AAGLTAAYAR                           73831000     73408000     74801000
      APEAEQVLSAAATFPIAQPATDVEAR           32425000     34162000     30033000
      AVQDNGESAFR                          15753000     14619000     16620000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        16317000     13280000     16185000
      GFTLAEVK                             35869000     32868000     36511000
      IAPRPLDLLRPVVR                       22705000     12422000     13232000
      IIVFPR                               66043000     62327000     56243000
      NQEIFDANVQR                          94201000     97850000     97828000
      TIGIAVDHR                                  NA     27610000     30487000
      VHFDQAGK                             21897000     25643000     24551000
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           255170000    239600000    232950000
      ANGTTVLVGMPAGAK                      24783000     21088000     22228000
      ATDGGAHGVINVSVSEAAIEASTR            115520000    138720000    170230000
      CCSDVFNQVVK                          56947000     59907000     55909000
      DIVGAVLK                             63901000     62136000     62045000
      EALDFFAR                            224370000    215930000    209770000
      EKDIVGAVLK                           47384000     48729000     44169000
      GVIFYESHGK                          128450000    136440000    150730000
      IGDYAGIK                            170190000    169400000    169190000
      LPLVGGHEGAGVVVGMGENVK               159280000    162300000    131700000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         22408000     20918000     39602000
      SIGGEVFIDFTK                         63157000     68817000     51520000
      SIPETQK                              15143000     11129000     10265000
      SISIVGSYVGNR                        347960000    336540000    320680000
      VLGIDGGEGK                           26500000     27786000     30555000
      VLGIDGGEGKEELFR                     376900000    366510000    342830000
      VVGLSTLPEIYEK                       591740000    509800000    563150000
      YSGVCHTDLHAWHGDWPLPVK               120300000    105060000    118290000
      ANGTVVLVGLPAGAK                       6666500      6713300      4377700
      CSSDVFNHVVK                           4862700      4156600      4399100
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        39762000     40097000     40341000
      DIPVPEPKPNEILINVK                    24035000     22053000     20636000
      EALDFFSR                              8182900      8111100           NA
      GVIFYENK                                   NA      7659100      5123000
      IQQGTDLAEVAPILCAGVTVYK                7828700      7746200      6872100
      IVGLSELPK                            24640000     21307000     22154000
      NMVSDIQEATK                           6084700      4691700      6858900
      VLGIDAGEEK                                 NA           NA           NA
                                       05000amol_R1 05000amol_R2 05000amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                            NA      2245100           NA
      AEWALR                                5855600      4951000      5088000
      DEGLHTDFACLLFAHLK                     5431700           NA           NA
      DIHDWNNR                                   NA           NA           NA
      ELETLREENR                                 NA      4476600           NA
      ESEFLFNAIHTIPEIGEK                   29365000     32819000     34116000
      GMMPGLTFSNELICR                      12050000     10131000      8760800
      IVTEAVEIEQR                          16551000     16858000     18016000
      LLVAFGNK                             10618000           NA     11095000
      LLVAFGNKK                                  NA      1405200           NA
      NKPDPAIVEK                           16009000     16631000     17552000
      TNFFEK                                6488500      7516200      6069900
      TVLFPIK                              11783000     12212000     10672000
      VENPFDFMENISLAGK                     19204000     16765000     18244000
      WIQDADALFGER                         19638000     21174000     19946000
      YFLDALPVALLGMNADLMNQYVEFVADR         23475000     23032000     23512000
      AANLGGVAVSGLEMAQNSQK                  9171700      5909800      9607500
      DAVWFGPPK                           163980000           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           5547700      4855100      4928900
      FLGFEQIFK                            41632000     43747000     38452000
      GANIASFVMVADAMLDQGDVF                54392000     50645000     53524000
      GCIISETGITSEQIHDIASAK                 7523500      7752000      6725200
      GGLCVDLK                             11188000           NA      8928500
      ICYAFMR                               9772700      9965400      8956800
      NSWEGVLTGK                           10818000     12374000     12567000
      SLEEIVDEYSTFSESK                      9056500      8346000      8665300
      VLPIVSVPER                           44170000     45348000     42537000
      VTISGSGNVAQYAALK                           NA      6520200      2968900
      VTWENDNGEQEVAQGYR                     3797200      3316100      3540200
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 36478000     30917000      8902700
      ALVAQGVK                             11876000     13251000     12316000
      FIAEGSNMGSTPEAIAVFETAR               41019000     39142000     37644000
      GANIASFIK                            19327000     16584000     18711000
      GCIISETGITSEQVADISSAK                 4437000      4159300      2672800
      HIGQDTDVPAGDIGVGGR                   33222000     31486000     27562000
      IMINCFNECIDYAK                       15101000     11838000     13346000
      ITWTSER                              17860000     21368000     18150000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        81966000     93082000     78985000
      SLEQIVNEYSTFSENK                     63205000     58515000     54849000
      STATGPSEAVWYGPPK                     44663000     48340000     49658000
      VDIALPCATQNEVSGEEAK                  53196000     49142000     54984000
      VIELGGTVVSLSDSK                      13868000     14514000     13164000
      VQYIAGARPWTHVQK                      13966000      4334500      9447700
      VTWENDKGEQEVAQGYR                    34863000     35762000     34957000
      AAGLTAAYAR                           56674000     59860000     56589000
      APEAEQVLSAAATFPIAQPATDVEAR           28597000     27560000     27790000
      AVQDNGESAFR                          14360000     14340000     12363000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        32396000     28015000     27295000
      GFTLAEVK                             28625000     30395000     28597000
      IAPRPLDLLRPVVR                       19446000     16928000      9539600
      IIVFPR                               55047000     54330000     55103000
      NQEIFDANVQR                          97314000     88930000     93185000
      TIGIAVDHR                            28376000     28133000     25086000
      VHFDQAGK                              9480100      9404800      9381100
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           165120000    167430000    161860000
      ANGTTVLVGMPAGAK                      22431000     18316000     20300000
      ATDGGAHGVINVSVSEAAIEASTR             83974000     89722000     96564000
      CCSDVFNQVVK                          49537000     50409000     49880000
      DIVGAVLK                                   NA     54320000           NA
      EALDFFAR                            198790000    201530000    192620000
      EKDIVGAVLK                           31788000     34168000     34761000
      GVIFYESHGK                           89159000     87110000     93689000
      IGDYAGIK                            139160000    145850000    131210000
      LPLVGGHEGAGVVVGMGENVK               126630000    129240000    128900000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         25320000     14734000     22950000
      SIGGEVFIDFTK                         66386000     65128000     59226000
      SIPETQK                              35061000     33558000     18176000
      SISIVGSYVGNR                        272760000    274860000    271630000
      VLGIDGGEGK                           19514000     20861000     20564000
      VLGIDGGEGKEELFR                     300540000    305140000    301640000
      VVGLSTLPEIYEK                       542550000    536350000    543290000
      YSGVCHTDLHAWHGDWPLPVK               123800000    136520000    148650000
      ANGTVVLVGLPAGAK                       4174600      4888700      4201100
      CSSDVFNHVVK                                NA      3599100      1324300
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        40341000     41004000     38382000
      DIPVPEPKPNEILINVK                    19922000     21247000     20546000
      EALDFFSR                              6808500      6789300      6528400
      GVIFYENK                                   NA           NA           NA
      IQQGTDLAEVAPILCAGVTVYK                7543700      9067000      8002600
      IVGLSELPK                            18601000     20152000     18739000
      NMVSDIQEATK                           5489200      5901200      7053100
      VLGIDAGEEK                                 NA           NA      4374200
                                       00500amol_R1 00500amol_R2 00500amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                      13581000      3857100      3408400
      AEWALR                                5981800      5561300      5153400
      DEGLHTDFACLLFAHLK                          NA           NA           NA
      DIHDWNNR                              2065700           NA           NA
      ELETLREENR                            4544500      3601600      6076700
      ESEFLFNAIHTIPEIGEK                   33732000     36075000     33870000
      GMMPGLTFSNELICR                      11000000      8619700      8180300
      IVTEAVEIEQR                          19916000     18134000     18131000
      LLVAFGNK                              7850300      9171900      9439800
      LLVAFGNKK                             2434900           NA           NA
      NKPDPAIVEK                           19665000     20910000     21316000
      TNFFEK                                7301400      5672100      6357000
      TVLFPIK                              10149000     10174000     12701000
      VENPFDFMENISLAGK                     18620000     15523000     17782000
      WIQDADALFGER                         17870000     15877000     15662000
      YFLDALPVALLGMNADLMNQYVEFVADR         28671000     26459000     24449000
      AANLGGVAVSGLEMAQNSQK                 11134000     10031000     10783000
      DAVWFGPPK                                  NA           NA    173600000
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           5049900      5412800      4896300
      FLGFEQIFK                            44682000     49538000     45829000
      GANIASFVMVADAMLDQGDVF                73304000     71656000     72863000
      GCIISETGITSEQIHDIASAK                 7035500      9470300      8603100
      GGLCVDLK                             11027000     10436000           NA
      ICYAFMR                              12064000     13340000     11924000
      NSWEGVLTGK                           14270000     12867000     15402000
      SLEEIVDEYSTFSESK                     10537000      8580900      8875100
      VLPIVSVPER                           42713000     45913000     42373000
      VTISGSGNVAQYAALK                           NA           NA       114880
      VTWENDNGEQEVAQGYR                     3339600      2656700      2754200
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 35660000      9801500     37708000
      ALVAQGVK                             13475000     14918000     14473000
      FIAEGSNMGSTPEAIAVFETAR               41312000     44045000     42909000
      GANIASFIK                            21582000     22035000     22899000
      GCIISETGITSEQVADISSAK                 4240700      3561200      3530800
      HIGQDTDVPAGDIGVGGR                   36389000     43352000     34656000
      IMINCFNECIDYAK                       13809000     14127000     14172000
      ITWTSER                              18777000     20111000     20650000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        98332000    100620000    103800000
      SLEQIVNEYSTFSENK                     57182000     52723000     52679000
      STATGPSEAVWYGPPK                     45712000     50203000     45932000
      VDIALPCATQNEVSGEEAK                  60393000     56938000     55713000
      VIELGGTVVSLSDSK                      15646000     16003000     13963000
      VQYIAGARPWTHVQK                       6546500      9333900           NA
      VTWENDKGEQEVAQGYR                    41604000     41930000     40345000
      AAGLTAAYAR                           67658000     68629000     68448000
      APEAEQVLSAAATFPIAQPATDVEAR           20302000     18780000     14549000
      AVQDNGESAFR                          13298000     14881000     15462000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        25085000     29607000     27163000
      GFTLAEVK                             32146000     35644000     35608000
      IAPRPLDLLRPVVR                       31518000     24931000     29313000
      IIVFPR                               50827000     52423000     50617000
      NQEIFDANVQR                          98888000     92385000     90198000
      TIGIAVDHR                            23862000     23404000     25101000
      VHFDQAGK                             12013000     11715000     11022000
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           183080000    170030000    152060000
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             86872000     96722000    100900000
      CCSDVFNQVVK                          50684000     56791000     53475000
      DIVGAVLK                             61411000     49703000     51248000
      EALDFFAR                            215740000    217940000    209560000
      EKDIVGAVLK                           38526000     38143000     34845000
      GVIFYESHGK                          101610000    102340000     98705000
      IGDYAGIK                            131910000    127840000    136180000
      LPLVGGHEGAGVVVGMGENVK               149250000    123430000    151460000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         40856000     31000000     37518000
      SIGGEVFIDFTK                         58339000     70890000     80179000
      SIPETQK                              25544000     21280000     20038000
      SISIVGSYVGNR                        288100000    290000000    278750000
      VLGIDGGEGK                           10734000      8886300     10019000
      VLGIDGGEGKEELFR                     341350000    367580000    364430000
      VVGLSTLPEIYEK                       551250000    604480000    555250000
      YSGVCHTDLHAWHGDWPLPVK               129790000     92341000    119170000
      ANGTVVLVGLPAGAK                       5176000      5871100      5028300
      CSSDVFNHVVK                           4628900           NA      1072600
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        37040000     34270000     35832000
      DIPVPEPKPNEILINVK                          NA     20424000     22525000
      EALDFFSR                              7007000      6903200     10314000
      GVIFYENK                              5495200      5192100      6189600
      IQQGTDLAEVAPILCAGVTVYK                9206400      8705200      8108000
      IVGLSELPK                            20179000     22906000     22137000
      NMVSDIQEATK                           6515100      6892600      4909200
      VLGIDAGEEK                                 NA           NA           NA
                                       00050amol_R1 00050amol_R2 00050amol_R3
      AAADALSDLEIK                         20886000           NA           NA
      AAADALSDLEIKDSK                       2214900           NA           NA
      AEWALR                                5747200      6105600      5694800
      DEGLHTDFACLLFAHLK                    13861000      8548300      6367100
      DIHDWNNR                                   NA      2247900      2202500
      ELETLREENR                            8499500      8460600      7902200
      ESEFLFNAIHTIPEIGEK                   38346000     36535000     29781000
      GMMPGLTFSNELICR                       6151900           NA           NA
      IVTEAVEIEQR                          30619000     25508000     23114000
      LLVAFGNK                              8127400      7157900      7904400
      LLVAFGNKK                             2362900      2802500      3551400
      NKPDPAIVEK                           22673000     18246000     22252000
      TNFFEK                                6226000      7341400      6135600
      TVLFPIK                              10629000     10796000     11721000
      VENPFDFMENISLAGK                     14715000     15365000     15861000
      WIQDADALFGER                         17249000     14802000     15931000
      YFLDALPVALLGMNADLMNQYVEFVADR         33055000     31808000     29785000
      AANLGGVAVSGLEMAQNSQK                       NA     11810000      8983300
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                           32328000     25658000     28783000
      FHPSVNLSILK                          14125000      8598000      8592100
      FLGFEQIFK                            65732000     61723000     58388000
      GANIASFVMVADAMLDQGDVF               132300000     99299000     92946000
      GCIISETGITSEQIHDIASAK                10779000     11435000      8390700
      GGLCVDLK                             12594000     11354000     12042000
      ICYAFMR                              11076000     10936000     10639000
      NSWEGVLTGK                           13001000     15028000     15520000
      SLEEIVDEYSTFSESK                     11166000      8548800      6445100
      VLPIVSVPER                           36438000     34005000     36401000
      VTISGSGNVAQYAALK                      4470300      3338400           NA
      VTWENDNGEQEVAQGYR                     3588500      3160500      3463000
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK      9561000      6181500      5222300
      AANLGGVAVSGLEMAQNSQR                 32679000     30286000     31639000
      ALVAQGVK                                   NA       756380           NA
      FIAEGSNMGSTPEAIAVFETAR               44255000     43687000     42772000
      GANIASFIK                            17833000     19129000     18912000
      GCIISETGITSEQVADISSAK                 8135900      5579800           NA
      HIGQDTDVPAGDIGVGGR                   48231000     39830000     36437000
      IMINCFNECIDYAK                       14038000     15355000     14142000
      ITWTSER                                    NA     20222000     18486000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR       101350000           NA    106000000
      SLEQIVNEYSTFSENK                     53812000     10437000     48862000
      STATGPSEAVWYGPPK                     46355000     41778000     41345000
      VDIALPCATQNEVSGEEAK                   9008400     48730000     49888000
      VIELGGTVVSLSDSK                      20895000     16526000     14645000
      VQYIAGARPWTHVQK                            NA           NA           NA
      VTWENDKGEQEVAQGYR                    65052000     54825000     53014000
      AAGLTAAYAR                           55470000     62138000     59012000
      APEAEQVLSAAATFPIAQPATDVEAR            5320200      3008100      3149500
      AVQDNGESAFR                          18130000     18466000     15488000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        33564000     25133000     32184000
      GFTLAEVK                             34750000     37621000     36671000
      IAPRPLDLLRPVVR                       21233000     10944000     10565000
      IIVFPR                                2242100      2204500     46295000
      NQEIFDANVQR                          81259000     85991000     84943000
      TIGIAVDHR                            35365000     32625000     27458000
      VHFDQAGK                             14072000     15904000     15879000
      VHFDQAGKK                             2837100      2097300      1809700
      ANELLINVK                           184930000    178880000    182820000
      ANGTTVLVGMPAGAK                      22679000     20585000     20851000
      ATDGGAHGVINVSVSEAAIEASTR            124550000    101040000     90585000
      CCSDVFNQVVK                          71812000     66607000     69073000
      DIVGAVLK                             60042000     61109000     58708000
      EALDFFAR                            227080000    210410000    212330000
      EKDIVGAVLK                           62427000     35459000     33921000
      GVIFYESHGK                          229190000    226950000    171090000
      IGDYAGIK                                   NA    129730000    128600000
      LPLVGGHEGAGVVVGMGENVK               128030000    110050000    136420000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         61716000     46075000     44490000
      SIGGEVFIDFTK                         89599000     44375000     65134000
      SIPETQK                              33610000     32138000     59548000
      SISIVGSYVGNR                        299540000    265510000    271690000
      VLGIDGGEGK                           11276000      9520100      9843400
      VLGIDGGEGKEELFR                     453940000    460560000    487960000
      VVGLSTLPEIYEK                       563450000    526950000    557180000
      YSGVCHTDLHAWHGDWPLPVK                84953000    124150000    106120000
      ANGTVVLVGLPAGAK                       6290300      5874900      5615000
      CSSDVFNHVVK                           6555400      4441900      5483200
      DIPVPKPKPNELLINVK                     8572000      6594800      4489600
      VVGLSSLPEIYEK                        34226000     31966000     35836000
      DIPVPEPKPNEILINVK                     4256400     28021000     21048000
      EALDFFSR                              7828000      6971000      6133400
      GVIFYENK                              7300500      7597400      6775400
      IQQGTDLAEVAPILCAGVTVYK                7743900      8589700      8480300
      IVGLSELPK                            16650000     17530000     17048000
      NMVSDIQEATK                           7612200      7759200      7525400
      VLGIDAGEEK                                 NA           NA           NA

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
                         sample
      12500amol_R1 12500amol_R1
      12500amol_R2 12500amol_R2
      12500amol_R3 12500amol_R3
      00125amol_R1 00125amol_R1
      00125amol_R2 00125amol_R2
      00125amol_R3 00125amol_R3
      25000amol_R1 25000amol_R1
      25000amol_R2 25000amol_R2
      25000amol_R3 25000amol_R3
      02500amol_R1 02500amol_R1
      02500amol_R2 02500amol_R2
      02500amol_R3 02500amol_R3
      00250amol_R1 00250amol_R1
      00250amol_R2 00250amol_R2
      00250amol_R3 00250amol_R3
      50000amol_R1 50000amol_R1
      50000amol_R2 50000amol_R2
      50000amol_R3 50000amol_R3
      05000amol_R1 05000amol_R1
      05000amol_R2 05000amol_R2
      05000amol_R3 05000amol_R3
      00500amol_R1 00500amol_R1
      00500amol_R2 00500amol_R2
      00500amol_R3 00500amol_R3
      00050amol_R1 00050amol_R1
      00050amol_R2 00050amol_R2
      00050amol_R3 00050amol_R3

---

    Code
      D2
    Output
      class: SummarizedExperiment 
      dim: 87 27 
      metadata(0):
      assays(1): intensities
      rownames(87): AAADALSDLEIK AAADALSDLEIKDSK ... NMVSDIQEATK VLGIDAGEEK
      rowData names(3): Sequence Proteins Score
      colnames(27): 12500amol_R1 12500amol_R2 ... 00050amol_R2 00050amol_R3
      colData names(1): sample

---

    Code
      SummarizedExperiment::assays(D2)$intensities
    Output
                                       12500amol_R1 12500amol_R2 12500amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                       8047500      7828500      5176800
      AEWALR                                5203200      4796100      5020600
      DEGLHTDFACLLFAHLK                     4653600      3672600       488180
      DIHDWNNR                              1983600      1854900           NA
      ELETLREENR                                 NA           NA           NA
      ESEFLFNAIHTIPEIGEK                   34234000     36632000     33268000
      GMMPGLTFSNELICR                       9125500      7441000      8216700
      IVTEAVEIEQR                          19044000     19817000     19294000
      LLVAFGNK                             11947000     12955000     11914000
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                           16866000     16956000     17424000
      TNFFEK                                7189600      7309700      5015900
      TVLFPIK                              14240000     13796000     12932000
      VENPFDFMENISLAGK                     17223000     18800000     19595000
      WIQDADALFGER                         18108000     19518000     18530000
      YFLDALPVALLGMNADLMNQYVEFVADR         26632000     30880000     28603000
      AANLGGVAVSGLEMAQNSQK                  9531300      3603100     10011000
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           4478600           NA      4808900
      FLGFEQIFK                            41513000     41062000     40552000
      GANIASFVMVADAMLDQGDVF                57811000     56199000     56801000
      GCIISETGITSEQIHDIASAK                 8618800      8428600      8425500
      GGLCVDLK                                   NA     10517000           NA
      ICYAFMR                               9595500      9692600      9753000
      NSWEGVLTGK                           12479000     12264000     14467000
      SLEEIVDEYSTFSESK                      8899400      8597900     10299000
      VLPIVSVPER                           41015000     37670000     38551000
      VTISGSGNVAQYAALK                      4953700       139030      5124400
      VTWENDNGEQEVAQGYR                     3666600           NA      3181200
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 29672000     32350000     27713000
      ALVAQGVK                             11333000     13517000     12532000
      FIAEGSNMGSTPEAIAVFETAR               24002000     41679000     43253000
      GANIASFIK                            19570000     16194000     17111000
      GCIISETGITSEQVADISSAK                 3409200      3793000      3773900
      HIGQDTDVPAGDIGVGGR                   28595000     32818000     33506000
      IMINCFNECIDYAK                       14511000     14191000     12607000
      ITWTSER                              20711000     18026000     16717000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        11102000     89102000     92008000
      SLEQIVNEYSTFSENK                     55353000     69945000     56380000
      STATGPSEAVWYGPPK                     50388000     40574000     39416000
      VDIALPCATQNEVSGEEAK                  49569000     46119000     49729000
      VIELGGTVVSLSDSK                      14031000     15112000     16164000
      VQYIAGARPWTHVQK                            NA      6834000      4173500
      VTWENDKGEQEVAQGYR                    34850000     33313000     37568000
      AAGLTAAYAR                           65040000     61449000     52189000
      APEAEQVLSAAATFPIAQPATDVEAR           29710000     30017000     28594000
      AVQDNGESAFR                          10871000     13576000     12546000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        28071000     28998000     28216000
      GFTLAEVK                             29152000     30119000     28740000
      IAPRPLDLLRPVVR                       17729000     17374000     19096000
      IIVFPR                               51639000     48895000     49167000
      NQEIFDANVQR                          87065000     88759000     83979000
      TIGIAVDHR                            25122000     27070000     23807000
      VHFDQAGK                              9541900      9842900     11953000
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           166710000    177430000    176800000
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR            112140000    112550000    118760000
      CCSDVFNQVVK                          51288000     50276000     47376000
      DIVGAVLK                                   NA     52654000           NA
      EALDFFAR                            191980000    190470000    185860000
      EKDIVGAVLK                           39383000     38553000     40919000
      GVIFYESHGK                           89749000     95686000    101170000
      IGDYAGIK                            143460000    130830000    137570000
      LPLVGGHEGAGVVVGMGENVK               121990000    123740000    142930000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         28729000     30246000     23600000
      SIGGEVFIDFTK                         72612000     66860000     56807000
      SIPETQK                              30419000     27045000     26483000
      SISIVGSYVGNR                        263390000    281730000    259550000
      VLGIDGGEGK                           21975000     21132000           NA
      VLGIDGGEGKEELFR                     294480000    305320000    314990000
      VVGLSTLPEIYEK                       538290000    508780000    520300000
      YSGVCHTDLHAWHGDWPLPVK               100260000     93705000    102610000
      ANGTVVLVGLPAGAK                       5215500      4563500      4349800
      CSSDVFNHVVK                                NA           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        40390000     36848000     37237000
      DIPVPEPKPNEILINVK                    18878000     21630000     17152000
      EALDFFSR                              7604100      5961700      6189500
      GVIFYENK                              5623000      5376300      5046900
      IQQGTDLAEVAPILCAGVTVYK                8343800      8187400      5859200
      IVGLSELPK                            17116000     16903000     18721000
      NMVSDIQEATK                           5862500      6922900      5567500
      VLGIDAGEEK                                 NA           NA           NA
                                       00125amol_R1 00125amol_R2 00125amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                      12535000     12650000     12914000
      AEWALR                                5886400      5541000      5128600
      DEGLHTDFACLLFAHLK                     2918700      4592400      4020300
      DIHDWNNR                              2016200           NA      2459200
      ELETLREENR                            5044000      4592700      6420100
      ESEFLFNAIHTIPEIGEK                   32775000     30116000     38596000
      GMMPGLTFSNELICR                       6410400      7156500      7707300
      IVTEAVEIEQR                          17074000     17950000     18123000
      LLVAFGNK                              6318600      7032200      7450800
      LLVAFGNKK                             1617900      1922000      1758300
      NKPDPAIVEK                           18798000     17543000     16623000
      TNFFEK                                5626500      4881700      4805400
      TVLFPIK                              11262000      9855900     10123000
      VENPFDFMENISLAGK                     18194000     16374000     14859000
      WIQDADALFGER                         16198000     15526000     16581000
      YFLDALPVALLGMNADLMNQYVEFVADR         24598000     21163000     20556000
      AANLGGVAVSGLEMAQNSQK                  9637200      8951900      9468500
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                           23975000     24214000     23049000
      FHPSVNLSILK                           4761100      4566800      4932100
      FLGFEQIFK                            48344000     42988000     45800000
      GANIASFVMVADAMLDQGDVF                78191000     78934000     77210000
      GCIISETGITSEQIHDIASAK                 8114600      7191500      8416900
      GGLCVDLK                             10639000           NA           NA
      ICYAFMR                              11386000     12635000     11527000
      NSWEGVLTGK                           12799000     16737000           NA
      SLEEIVDEYSTFSESK                      9072200      9319500      8774500
      VLPIVSVPER                           32628000     36498000     31489000
      VTISGSGNVAQYAALK                      6602200       107990           NA
      VTWENDNGEQEVAQGYR                          NA      2820200      2282200
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA      1790400
      AANLGGVAVSGLEMAQNSQR                  9409000     33196000     33261000
      ALVAQGVK                             14090000     13305000     13422000
      FIAEGSNMGSTPEAIAVFETAR               39925000     39763000     42135000
      GANIASFIK                            18180000     17768000     18940000
      GCIISETGITSEQVADISSAK                 2851800      2931900      2648400
      HIGQDTDVPAGDIGVGGR                   33815000     33476000     32333000
      IMINCFNECIDYAK                       12242000     12196000     12046000
      ITWTSER                              16961000     16113000       515960
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR         9245600    105750000     96463000
      SLEQIVNEYSTFSENK                     45184000     48337000     43716000
      STATGPSEAVWYGPPK                     37407000     40732000     45507000
      VDIALPCATQNEVSGEEAK                  45972000     49351000     48412000
      VIELGGTVVSLSDSK                      12814000     12436000     12776000
      VQYIAGARPWTHVQK                       3603500           NA      3652800
      VTWENDKGEQEVAQGYR                    42150000     38839000     43270000
      AAGLTAAYAR                           57292000     56821000     59909000
      APEAEQVLSAAATFPIAQPATDVEAR           19276000     23141000     20579000
      AVQDNGESAFR                          14132000     12807000     14194000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        39555000     40572000     40075000
      GFTLAEVK                             31231000     31843000     31087000
      IAPRPLDLLRPVVR                       22593000     26704000     26954000
      IIVFPR                               40028000     41502000      3146600
      NQEIFDANVQR                          77218000     74951000     76925000
      TIGIAVDHR                            19422000     18351000     16928000
      VHFDQAGK                                   NA           NA     13287000
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           156250000    162970000    160710000
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             78673000     59114000     77899000
      CCSDVFNQVVK                          54965000     55758000     52645000
      DIVGAVLK                                   NA     54211000     49563000
      EALDFFAR                            187870000    188160000    190960000
      EKDIVGAVLK                           30932000     31701000     33596000
      GVIFYESHGK                          111700000    108790000    104130000
      IGDYAGIK                            128550000    124470000    132530000
      LPLVGGHEGAGVVVGMGENVK               107330000     78272000    142440000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         22033000     27413000     31741000
      SIGGEVFIDFTK                         77111000     58227000     73145000
      SIPETQK                              26348000     44054000     46080000
      SISIVGSYVGNR                        239120000    256900000    261340000
      VLGIDGGEGK                            9609600     11221000     10503000
      VLGIDGGEGKEELFR                     383050000    355750000    310680000
      VVGLSTLPEIYEK                       500750000    499990000    510720000
      YSGVCHTDLHAWHGDWPLPVK                86925000     97304000    111570000
      ANGTVVLVGLPAGAK                       4540500      5673900      4842000
      CSSDVFNHVVK                           3969400      4748900      4075800
      DIPVPKPKPNELLINVK                          NA      1188800           NA
      VVGLSSLPEIYEK                        35070000     34282000     34474000
      DIPVPEPKPNEILINVK                    17538000     20456000     20083000
      EALDFFSR                              6564500      6941500      6239000
      GVIFYENK                                   NA      5730800      5494500
      IQQGTDLAEVAPILCAGVTVYK                7604500      7188900      7865100
      IVGLSELPK                            18118000     18386000     20462000
      NMVSDIQEATK                           3543100      4131200           NA
      VLGIDAGEEK                                 NA           NA           NA
                                       25000amol_R1 25000amol_R2 25000amol_R3
      AAADALSDLEIK                               NA           NA     31031000
      AAADALSDLEIKDSK                       7401800      8947400      5334500
      AEWALR                                4750000      5040800      5368400
      DEGLHTDFACLLFAHLK                     5027600      3054200      5661000
      DIHDWNNR                              2537200       594660      1966900
      ELETLREENR                                 NA       991370           NA
      ESEFLFNAIHTIPEIGEK                   32349000     27488000     34607000
      GMMPGLTFSNELICR                       7887300      8045400      7319400
      IVTEAVEIEQR                          18931000     20324000     22839000
      LLVAFGNK                             10795000     13952000     11907000
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                           15232000     15531000     19733000
      TNFFEK                                5658100      6115700      7342800
      TVLFPIK                              12106000     13714000     13279000
      VENPFDFMENISLAGK                     17741000     15629000     14853000
      WIQDADALFGER                         18736000     17277000     17286000
      YFLDALPVALLGMNADLMNQYVEFVADR         30426000     21571000     22501000
      AANLGGVAVSGLEMAQNSQK                  9042000      9672700      8450500
      DAVWFGPPK                           164440000           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           5824900      6678000      6571000
      FLGFEQIFK                            36591000     36452000     33500000
      GANIASFVMVADAMLDQGDVF                52637000     40938000     43741000
      GCIISETGITSEQIHDIASAK                 7610000      6353500      8576800
      GGLCVDLK                             10081000     10871000     11769000
      ICYAFMR                               8519100      6256900      8973500
      NSWEGVLTGK                           13811000     14103000     13184000
      SLEEIVDEYSTFSESK                      9935600      6729600      8751600
      VLPIVSVPER                           38915000     48419000     40096000
      VTISGSGNVAQYAALK                      1857800      4924100       955960
      VTWENDNGEQEVAQGYR                     2482000           NA      1882300
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 30710000     30190000     27414000
      ALVAQGVK                             12354000     16352000     17099000
      FIAEGSNMGSTPEAIAVFETAR               38018000     40932000     36395000
      GANIASFIK                            17563000     16593000           NA
      GCIISETGITSEQVADISSAK                 2982100      2920400      3289600
      HIGQDTDVPAGDIGVGGR                   29277000     32548000     35047000
      IMINCFNECIDYAK                       11376000     11134000     11175000
      ITWTSER                              17454000     20639000     20497000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        96838000     66771000     82458000
      SLEQIVNEYSTFSENK                     54289000     48529000     55413000
      STATGPSEAVWYGPPK                     41070000     48208000     40245000
      VDIALPCATQNEVSGEEAK                  48984000     51460000     48491000
      VIELGGTVVSLSDSK                      13786000     15033000     14579000
      VQYIAGARPWTHVQK                       3555600      6063400      2822500
      VTWENDKGEQEVAQGYR                    33313000     29170000     30885000
      AAGLTAAYAR                           56708000     64624000     67064000
      APEAEQVLSAAATFPIAQPATDVEAR           31078000     28494000     32203000
      AVQDNGESAFR                          11641000     12948000     14337000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        18878000     22654000     23250000
      GFTLAEVK                             29366000     27287000     27300000
      IAPRPLDLLRPVVR                       16956000     13445000     18669000
      IIVFPR                               50223000     57016000     50919000
      NQEIFDANVQR                          79339000     85649000     85974000
      TIGIAVDHR                            24726000     24372000     28093000
      VHFDQAGK                             11753000     10328000     12035000
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           175030000    179560000    171910000
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             81553000     84460000    109120000
      CCSDVFNQVVK                          45459000     48760000     47025000
      DIVGAVLK                             47512000     54924000     49717000
      EALDFFAR                            187090000    202670000    172850000
      EKDIVGAVLK                           37427000     30443000     31028000
      GVIFYESHGK                           92680000     88425000     94108000
      IGDYAGIK                            140860000    164480000    163970000
      LPLVGGHEGAGVVVGMGENVK               121090000     91603000    112810000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         27627000     41993000     31996000
      SIGGEVFIDFTK                         89950000     54981000     48927000
      SIPETQK                              23624000     20399000     14067000
      SISIVGSYVGNR                        268050000    271070000    266050000
      VLGIDGGEGK                           20146000     27612000     29202000
      VLGIDGGEGKEELFR                     298300000    254790000    265170000
      VVGLSTLPEIYEK                       492600000    515000000    484600000
      YSGVCHTDLHAWHGDWPLPVK               109910000     83720000    104270000
      ANGTVVLVGLPAGAK                       4344900      4760500      4648500
      CSSDVFNHVVK                           4227700           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        39022000     40275000     37821000
      DIPVPEPKPNEILINVK                    21322000     22581000     22346000
      EALDFFSR                              6284000      6514300      5086900
      GVIFYENK                              4564900      4010700      6011100
      IQQGTDLAEVAPILCAGVTVYK                8015500      6748800      6188500
      IVGLSELPK                            18717000     16968000     16518000
      NMVSDIQEATK                           5690000      6231300      6148900
      VLGIDAGEEK                                 NA           NA           NA
                                       02500amol_R1 02500amol_R2 02500amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                       6136200      7674900      8591400
      AEWALR                                5957600      4786100      4816700
      DEGLHTDFACLLFAHLK                     5267400      5164700      3786000
      DIHDWNNR                               670810           NA           NA
      ELETLREENR                                 NA           NA      4589300
      ESEFLFNAIHTIPEIGEK                   36682000     26819000     35535000
      GMMPGLTFSNELICR                       8917400      9855100      9280200
      IVTEAVEIEQR                          18746000     18706000     18415000
      LLVAFGNK                             11477000     10808000      9501100
      LLVAFGNKK                              812840           NA           NA
      NKPDPAIVEK                           15359000     15020000     18017000
      TNFFEK                                7093800      5883800      6104500
      TVLFPIK                              10853000     11873000     12222000
      VENPFDFMENISLAGK                     17804000     16232000     16165000
      WIQDADALFGER                         20448000     22291000     19101000
      YFLDALPVALLGMNADLMNQYVEFVADR         26373000     24192000     24640000
      AANLGGVAVSGLEMAQNSQK                 11468000     10733000      8357800
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           5795700      4869500      7744400
      FLGFEQIFK                            42337000     43369000     39303000
      GANIASFVMVADAMLDQGDVF                53492000     57457000     57307000
      GCIISETGITSEQIHDIASAK                 6487200      8265500      8534100
      GGLCVDLK                             10390000      9708100     11243000
      ICYAFMR                               9011900      9064700      9718700
      NSWEGVLTGK                           13538000     13777000     11775000
      SLEEIVDEYSTFSESK                      8935700     11227000      8276900
      VLPIVSVPER                           43658000     40543000     42875000
      VTISGSGNVAQYAALK                      7304200       179200           NA
      VTWENDNGEQEVAQGYR                     3591000      2945800      3175700
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 33938000     34882000     32829000
      ALVAQGVK                                   NA     12826000     13388000
      FIAEGSNMGSTPEAIAVFETAR               38443000     38338000     36175000
      GANIASFIK                            19074000     18107000     18392000
      GCIISETGITSEQVADISSAK                 3964200      3266500      3577100
      HIGQDTDVPAGDIGVGGR                   36245000     29174000     32116000
      IMINCFNECIDYAK                       13813000     14105000     13601000
      ITWTSER                                    NA     17431000     16695000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        93842000     97643000     84389000
      SLEQIVNEYSTFSENK                     55704000     53728000     58928000
      STATGPSEAVWYGPPK                     51951000     50880000     44896000
      VDIALPCATQNEVSGEEAK                  56650000     57853000     54369000
      VIELGGTVVSLSDSK                      14685000     14603000     13641000
      VQYIAGARPWTHVQK                       8729600           NA      6308900
      VTWENDKGEQEVAQGYR                    36113000     32947000     33843000
      AAGLTAAYAR                           57938000     56413000     56256000
      APEAEQVLSAAATFPIAQPATDVEAR           28101000     28751000     28906000
      AVQDNGESAFR                          14091000     13395000     13574000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        29886000     23200000     21104000
      GFTLAEVK                             31813000     32480000     30047000
      IAPRPLDLLRPVVR                        9304100     20256000     10159000
      IIVFPR                               62763000     57768000     51498000
      NQEIFDANVQR                          90602000     88374000     87222000
      TIGIAVDHR                            23628000     23605000     23287000
      VHFDQAGK                              9722900           NA     10828000
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           162450000    148300000    153900000
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             75336000     93820000     89931000
      CCSDVFNQVVK                          55132000     52541000     53446000
      DIVGAVLK                             50133000           NA     49281000
      EALDFFAR                            204860000    200500000    196980000
      EKDIVGAVLK                           31296000     31366000     36342000
      GVIFYESHGK                           89749000     87816000     91649000
      IGDYAGIK                            144240000    130430000    142680000
      LPLVGGHEGAGVVVGMGENVK               106610000    116970000    153450000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         39722000     21790000     28665000
      SIGGEVFIDFTK                         55031000     72442000     73249000
      SIPETQK                              42448000     21952000     23849000
      SISIVGSYVGNR                        275740000    266050000    256290000
      VLGIDGGEGK                           20516000     19955000     21155000
      VLGIDGGEGKEELFR                     292470000    310500000    331300000
      VVGLSTLPEIYEK                       554590000    524410000    503030000
      YSGVCHTDLHAWHGDWPLPVK               154610000    153570000    142580000
      ANGTVVLVGLPAGAK                       5258200      4035000      5106900
      CSSDVFNHVVK                                NA           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        38612000     39130000     34649000
      DIPVPEPKPNEILINVK                    23169000     24205000     22123000
      EALDFFSR                              7582400      6551000      8855600
      GVIFYENK                              5259400      5296400      5683000
      IQQGTDLAEVAPILCAGVTVYK                7065700      8628700      6986100
      IVGLSELPK                            20302000     19966000     19687000
      NMVSDIQEATK                           6636100      6956800      6133800
      VLGIDAGEEK                                 NA           NA           NA
                                       00250amol_R1 00250amol_R2 00250amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                      10380000     13746000     12170000
      AEWALR                                5012800      6114900      4947400
      DEGLHTDFACLLFAHLK                          NA           NA           NA
      DIHDWNNR                                   NA      1825600      2114800
      ELETLREENR                            4025500      4004000      3990500
      ESEFLFNAIHTIPEIGEK                   27081000     24408000     36055000
      GMMPGLTFSNELICR                       7410700           NA      7451800
      IVTEAVEIEQR                          17715000     17130000     18822000
      LLVAFGNK                              6810300      6424100      7090500
      LLVAFGNKK                                  NA      1740700           NA
      NKPDPAIVEK                           16117000     19753000     18334000
      TNFFEK                                6529300      5708000      5312500
      TVLFPIK                               9138900      9832600     11313000
      VENPFDFMENISLAGK                     17183000     17946000     17716000
      WIQDADALFGER                         14187000     14929000     16711000
      YFLDALPVALLGMNADLMNQYVEFVADR         23393000     23187000     23672000
      AANLGGVAVSGLEMAQNSQK                 10497000     10396000      5805700
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA     20789000     21193000
      FHPSVNLSILK                           4717400      6559500      5245400
      FLGFEQIFK                            45668000     45045000     43176000
      GANIASFVMVADAMLDQGDVF                73167000     83465000     75361000
      GCIISETGITSEQIHDIASAK                 8476100      8914100      6826000
      GGLCVDLK                             11181000     11518000     10404000
      ICYAFMR                              10124000     10851000     11857000
      NSWEGVLTGK                           11992000     13991000     13850000
      SLEEIVDEYSTFSESK                      8527400      8284100      7982700
      VLPIVSVPER                           32525000     33201000     38620000
      VTISGSGNVAQYAALK                      3240200           NA       156250
      VTWENDNGEQEVAQGYR                     2868800      2686600      3136100
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 35789000     34526000     38517000
      ALVAQGVK                             13542000     13044000     11838000
      FIAEGSNMGSTPEAIAVFETAR               41508000     27338000     41691000
      GANIASFIK                            17281000     18179000     18688000
      GCIISETGITSEQVADISSAK                 4058200      3657600      3434100
      HIGQDTDVPAGDIGVGGR                   32401000     36486000     30704000
      IMINCFNECIDYAK                       14150000     13363000     12834000
      ITWTSER                              17054000     17019000     16202000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR       101740000     97755000     95275000
      SLEQIVNEYSTFSENK                     52612000     53220000     58148000
      STATGPSEAVWYGPPK                     45311000     44881000     37105000
      VDIALPCATQNEVSGEEAK                  47434000     42819000     49205000
      VIELGGTVVSLSDSK                      13370000     13250000     12850000
      VQYIAGARPWTHVQK                      10381000     10087000      4458900
      VTWENDKGEQEVAQGYR                    35518000     41844000     38159000
      AAGLTAAYAR                           60960000     57043000     57159000
      APEAEQVLSAAATFPIAQPATDVEAR           19621000     18375000     14559000
      AVQDNGESAFR                            269700     13867000     11392000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        26657000     26955000     28359000
      GFTLAEVK                             29878000     34205000     30478000
      IAPRPLDLLRPVVR                       12654000     24344000     23507000
      IIVFPR                               42278000     47828000     42543000
      NQEIFDANVQR                          87783000     81450000     75827000
      TIGIAVDHR                            20533000     23004000     20200000
      VHFDQAGK                             10571000     10919000     10085000
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           159080000    167590000    153580000
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             75309000     76271000     63890000
      CCSDVFNQVVK                          50214000     56420000     54110000
      DIVGAVLK                             51955000     54844000           NA
      EALDFFAR                            184980000    181450000    197260000
      EKDIVGAVLK                           33061000     36486000     36213000
      GVIFYESHGK                           96308000     90249000     97108000
      IGDYAGIK                            130140000    128400000           NA
      LPLVGGHEGAGVVVGMGENVK               115130000    120610000     82944000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         26111000     30362000     28730000
      SIGGEVFIDFTK                         61958000     60162000     51938000
      SIPETQK                              27172000     23612000     24499000
      SISIVGSYVGNR                        237620000    249250000    250840000
      VLGIDGGEGK                            9815400     10868000     10356000
      VLGIDGGEGKEELFR                     306510000    323320000    334480000
      VVGLSTLPEIYEK                       516470000    485840000    510680000
      YSGVCHTDLHAWHGDWPLPVK                95987000    111420000    106180000
      ANGTVVLVGLPAGAK                       4616100      4648000      5237800
      CSSDVFNHVVK                                NA      3847200      4005300
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        35915000     33743000     33084000
      DIPVPEPKPNEILINVK                      740560     19561000     20542000
      EALDFFSR                              7199100           NA      6770800
      GVIFYENK                                   NA      5594700      5905600
      IQQGTDLAEVAPILCAGVTVYK                9389500      6423700      8971900
      IVGLSELPK                            17391000     18277000     19449000
      NMVSDIQEATK                           5782900      5473300      6290900
      VLGIDAGEEK                                 NA           NA           NA
                                       50000amol_R1 50000amol_R2 50000amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                       3703600      2446900      9229300
      AEWALR                                4289200      4344300      4488100
      DEGLHTDFACLLFAHLK                          NA      4253200      2208400
      DIHDWNNR                              2370700      1081100           NA
      ELETLREENR                            4575300      4366000      5064500
      ESEFLFNAIHTIPEIGEK                   33140000     31659000     27406000
      GMMPGLTFSNELICR                       7754600      6896300      5400100
      IVTEAVEIEQR                          21817000     21914000     20292000
      LLVAFGNK                             13680000     14222000     12507000
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                           18843000     20608000     20069000
      TNFFEK                                4697500      5281600      4999000
      TVLFPIK                              13729000     13299000     12643000
      VENPFDFMENISLAGK                     11665000     11872000     11978000
      WIQDADALFGER                         13153000     11702000     14635000
      YFLDALPVALLGMNADLMNQYVEFVADR         22177000     21612000     21594000
      AANLGGVAVSGLEMAQNSQK                  8643300      8194800      7711700
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                          11080000      5334300      6894900
      FLGFEQIFK                            34041000     36167000     34068000
      GANIASFVMVADAMLDQGDVF                37545000     37532000     39235000
      GCIISETGITSEQIHDIASAK                 6588900      6994200      8345100
      GGLCVDLK                             10352000     11960000     10222000
      ICYAFMR                               8104300      8011600      7782900
      NSWEGVLTGK                           13947000     13245000     12318000
      SLEEIVDEYSTFSESK                      4872000      6396000      4971500
      VLPIVSVPER                           39256000     37835000     33234000
      VTISGSGNVAQYAALK                      5208600      4380800           NA
      VTWENDNGEQEVAQGYR                          NA      2258800           NA
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 27270000     20562000     21709000
      ALVAQGVK                             13872000     13282000     12761000
      FIAEGSNMGSTPEAIAVFETAR               17129000     26417000     29255000
      GANIASFIK                            17478000     17191000     16950000
      GCIISETGITSEQVADISSAK                 3230700      3464700      4252600
      HIGQDTDVPAGDIGVGGR                   38034000     40391000     32215000
      IMINCFNECIDYAK                       13534000     10107000      9706900
      ITWTSER                              14538000     15590000           NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        66737000      9158900      8643600
      SLEQIVNEYSTFSENK                     47860000     42581000     48531000
      STATGPSEAVWYGPPK                     36156000           NA     35781000
      VDIALPCATQNEVSGEEAK                  47847000     45558000     39866000
      VIELGGTVVSLSDSK                      12599000     11492000     11850000
      VQYIAGARPWTHVQK                       3681500      3936800           NA
      VTWENDKGEQEVAQGYR                    30270000     34966000     35124000
      AAGLTAAYAR                           57364000     58738000     59985000
      APEAEQVLSAAATFPIAQPATDVEAR           25193000     27335000     24084000
      AVQDNGESAFR                          12240000     11697000     13328000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        12678000     10626000     12979000
      GFTLAEVK                             27870000     26299000     29279000
      IAPRPLDLLRPVVR                       17641000      9939700     10611000
      IIVFPR                               51313000     49871000     45102000
      NQEIFDANVQR                          73191000     78296000     78450000
      TIGIAVDHR                                  NA     22093000     24448000
      VHFDQAGK                             17013000     20518000     19688000
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           198260000    191720000    186800000
      ANGTTVLVGMPAGAK                      19256000     16874000     17825000
      ATDGGAHGVINVSVSEAAIEASTR             89756000    111000000    136510000
      CCSDVFNQVVK                          44246000     47935000     44835000
      DIVGAVLK                             49650000     49719000     49755000
      EALDFFAR                            174330000    172780000    168220000
      EKDIVGAVLK                           36816000     38991000     35420000
      GVIFYESHGK                           99804000    109170000    120870000
      IGDYAGIK                            132230000    135550000    135680000
      LPLVGGHEGAGVVVGMGENVK               123760000    129860000    105610000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         17410000     16738000     31758000
      SIGGEVFIDFTK                         49071000     55065000     41315000
      SIPETQK                              11766000      8904900      8231800
      SISIVGSYVGNR                        270350000    269280000    257160000
      VLGIDGGEGK                           20590000     22233000     24503000
      VLGIDGGEGKEELFR                     292840000    293270000    274920000
      VVGLSTLPEIYEK                       459770000    407920000    451610000
      YSGVCHTDLHAWHGDWPLPVK                93467000     84064000     94859000
      ANGTVVLVGLPAGAK                       5179700      5371700      3510600
      CSSDVFNHVVK                           3778200      3326000      3527800
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        30894000     32084000     32351000
      DIPVPEPKPNEILINVK                    18674000     17646000     16548000
      EALDFFSR                              6357900      6490200           NA
      GVIFYENK                                   NA      6128500      4108200
      IQQGTDLAEVAPILCAGVTVYK                6082700      6198200      5510900
      IVGLSELPK                            19145000     17049000     17766000
      NMVSDIQEATK                           4727600      3754100      5500300
      VLGIDAGEEK                                 NA           NA           NA
                                       05000amol_R1 05000amol_R2 05000amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                            NA      2319400           NA
      AEWALR                                6002800      5114900      5323300
      DEGLHTDFACLLFAHLK                     5568200           NA           NA
      DIHDWNNR                                   NA           NA           NA
      ELETLREENR                                 NA      4624700           NA
      ESEFLFNAIHTIPEIGEK                   30103000     33905000     35694000
      GMMPGLTFSNELICR                      12353000     10467000      9166100
      IVTEAVEIEQR                          16967000     17416000     18849000
      LLVAFGNK                             10885000           NA     11609000
      LLVAFGNKK                                  NA      1451700           NA
      NKPDPAIVEK                           16412000     17181000     18364000
      TNFFEK                                6651600      7765000      6350700
      TVLFPIK                              12079000     12616000     11166000
      VENPFDFMENISLAGK                     19687000     17320000     19087000
      WIQDADALFGER                         20132000     21875000     20868000
      YFLDALPVALLGMNADLMNQYVEFVADR         24065000     23794000     24599000
      AANLGGVAVSGLEMAQNSQK                  9402200      6105500     10052000
      DAVWFGPPK                           168100000           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           5687100      5015800      5156900
      FLGFEQIFK                            42678000     45195000     40231000
      GANIASFVMVADAMLDQGDVF                55759000     52321000     56000000
      GCIISETGITSEQIHDIASAK                 7712500      8008600      7036300
      GGLCVDLK                             11469000           NA      9341500
      ICYAFMR                              10018000     10295000      9371100
      NSWEGVLTGK                           11090000     12783000     13149000
      SLEEIVDEYSTFSESK                      9284100      8622200      9066200
      VLPIVSVPER                           45280000     46849000     44505000
      VTISGSGNVAQYAALK                           NA      6736000      3106300
      VTWENDNGEQEVAQGYR                     3892600      3425900      3704000
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 37395000     31941000      9314500
      ALVAQGVK                             12174000     13689000     12886000
      FIAEGSNMGSTPEAIAVFETAR               42050000     40438000     39385000
      GANIASFIK                            19813000     17133000     19576000
      GCIISETGITSEQVADISSAK                 4548600      4296900      2796400
      HIGQDTDVPAGDIGVGGR                   34057000     32529000     28836000
      IMINCFNECIDYAK                       15481000     12230000     13964000
      ITWTSER                              18309000     22076000     18989000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        84026000     96163000     82639000
      SLEQIVNEYSTFSENK                     64794000     60452000     57386000
      STATGPSEAVWYGPPK                     45786000     49940000     51955000
      VDIALPCATQNEVSGEEAK                  54533000     50768000     57527000
      VIELGGTVVSLSDSK                      14217000     14995000     13773000
      VQYIAGARPWTHVQK                      14317000      4478000      9884700
      VTWENDKGEQEVAQGYR                    35739000     36945000     36574000
      AAGLTAAYAR                           58098000     61842000     59206000
      APEAEQVLSAAATFPIAQPATDVEAR           29316000     28472000     29076000
      AVQDNGESAFR                          14721000     14815000     12935000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        33210000     28942000     28558000
      GFTLAEVK                             29344000     31401000     29920000
      IAPRPLDLLRPVVR                       19935000     17489000      9980900
      IIVFPR                               56431000     56129000     57652000
      NQEIFDANVQR                          99760000     91874000     97495000
      TIGIAVDHR                            29090000     29064000     26246000
      VHFDQAGK                              9718300      9716100      9815000
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           169270000    172980000    169350000
      ANGTTVLVGMPAGAK                      22994000     18922000     21239000
      ATDGGAHGVINVSVSEAAIEASTR             86084000     92691000    101030000
      CCSDVFNQVVK                          50782000     52078000     52187000
      DIVGAVLK                                   NA     56118000           NA
      EALDFFAR                            203790000    208200000    201530000
      EKDIVGAVLK                           32587000     35299000     36369000
      GVIFYESHGK                           91400000     89993000     98023000
      IGDYAGIK                            142660000    150680000    137280000
      LPLVGGHEGAGVVVGMGENVK               129820000    133520000    134870000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         25956000     15221000     24012000
      SIGGEVFIDFTK                         68054000     67284000     61965000
      SIPETQK                              35942000     34669000     19017000
      SISIVGSYVGNR                        279620000    283950000    284190000
      VLGIDGGEGK                           20004000     21552000     21515000
      VLGIDGGEGKEELFR                     308090000    315240000    315590000
      VVGLSTLPEIYEK                       556180000    554100000    568420000
      YSGVCHTDLHAWHGDWPLPVK               126910000    141040000    155530000
      ANGTVVLVGLPAGAK                       4279600      5050500      4395400
      CSSDVFNHVVK                                NA      3718200      1385500
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        41354000     42361000     40157000
      DIPVPEPKPNEILINVK                    20423000     21950000     21496000
      EALDFFSR                              6979600      7014000      6830400
      GVIFYENK                                   NA           NA           NA
      IQQGTDLAEVAPILCAGVTVYK                7733200      9367100      8372800
      IVGLSELPK                            19068000     20819000     19606000
      NMVSDIQEATK                           5627200      6096600      7379400
      VLGIDAGEEK                                 NA           NA      4576600
                                       00500amol_R1 00500amol_R2 00500amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                      12430000      3551700      3185400
      AEWALR                                5474900      5120900      4816200
      DEGLHTDFACLLFAHLK                          NA           NA           NA
      DIHDWNNR                              1890600           NA           NA
      ELETLREENR                            4159400      3316400      5679100
      ESEFLFNAIHTIPEIGEK                   30873000     33219000     31653000
      GMMPGLTFSNELICR                      10068000      7937100      7645000
      IVTEAVEIEQR                          18228000     16698000     16945000
      LLVAFGNK                              7185000      8445600      8822100
      LLVAFGNKK                             2228600           NA           NA
      NKPDPAIVEK                           17999000     19254000     19921000
      TNFFEK                                6682700      5222900      5941000
      TVLFPIK                               9288900      9368000     11870000
      VENPFDFMENISLAGK                     17042000     14294000     16618000
      WIQDADALFGER                         16356000     14620000     14637000
      YFLDALPVALLGMNADLMNQYVEFVADR         26241000     24364000     22849000
      AANLGGVAVSGLEMAQNSQK                 10190000      9236600     10077000
      DAVWFGPPK                                  NA           NA    162240000
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           4621900      4984200      4575900
      FLGFEQIFK                            40895000     45615000     42830000
      GANIASFVMVADAMLDQGDVF                67092000     65982000     68095000
      GCIISETGITSEQIHDIASAK                 6439300      8720400      8040200
      GGLCVDLK                             10093000      9609600           NA
      ICYAFMR                              11042000     12284000     11144000
      NSWEGVLTGK                           13061000     11848000     14395000
      SLEEIVDEYSTFSESK                      9644500      7901500      8294400
      VLPIVSVPER                           39093000     42277000     39600000
      VTISGSGNVAQYAALK                           NA           NA       107360
      VTWENDNGEQEVAQGYR                     3056600      2446300      2574000
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 32638000      9025300     35240000
      ALVAQGVK                             12333000     13737000     13526000
      FIAEGSNMGSTPEAIAVFETAR               37811000     40558000     40101000
      GANIASFIK                            19753000     20290000     21401000
      GCIISETGITSEQVADISSAK                 3881300      3279200      3299700
      HIGQDTDVPAGDIGVGGR                   33306000     39919000     32388000
      IMINCFNECIDYAK                       12639000     13008000     13244000
      ITWTSER                              17186000     18518000     19299000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        89999000     92656000     97006000
      SLEQIVNEYSTFSENK                     52336000     48548000     49232000
      STATGPSEAVWYGPPK                     41838000     46227000     42926000
      VDIALPCATQNEVSGEEAK                  55275000     52430000     52067000
      VIELGGTVVSLSDSK                      14320000     14735000     13049000
      VQYIAGARPWTHVQK                       5991700      8594800           NA
      VTWENDKGEQEVAQGYR                    38079000     38609000     37705000
      AAGLTAAYAR                           61924000     63195000     63969000
      APEAEQVLSAAATFPIAQPATDVEAR           18581000     17293000     13597000
      AVQDNGESAFR                          12171000     13703000     14450000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        22959000     27263000     25386000
      GFTLAEVK                             29422000     32822000     33278000
      IAPRPLDLLRPVVR                       28847000     22957000     27395000
      IIVFPR                               46519000     48272000     47305000
      NQEIFDANVQR                          90508000     85070000     84295000
      TIGIAVDHR                            21839000     21551000     23459000
      VHFDQAGK                             10995000     10787000     10301000
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           167560000    156560000    142110000
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             79510000     89063000     94301000
      CCSDVFNQVVK                          46389000     52294000     49976000
      DIVGAVLK                             56206000     45768000     47895000
      EALDFFAR                            197460000    200680000    195850000
      EKDIVGAVLK                           35261000     35123000     32565000
      GVIFYESHGK                           93002000     94236000     92246000
      IGDYAGIK                            120730000    117720000    127270000
      LPLVGGHEGAGVVVGMGENVK               136600000    113660000    141550000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         37393000     28546000     35063000
      SIGGEVFIDFTK                         53395000     65277000     74933000
      SIPETQK                              23379000     19595000     18727000
      SISIVGSYVGNR                        263690000    267040000    260500000
      VLGIDGGEGK                            9824200      8182700      9363600
      VLGIDGGEGKEELFR                     312430000    338470000    340580000
      VVGLSTLPEIYEK                       504530000    556620000    518920000
      YSGVCHTDLHAWHGDWPLPVK               118790000     85029000    111370000
      ANGTVVLVGLPAGAK                       4737300      5406200      4699300
      CSSDVFNHVVK                           4236600           NA      1002400
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        33901000     31556000     33487000
      DIPVPEPKPNEILINVK                          NA     18806000     21051000
      EALDFFSR                              6413200      6356600      9639200
      GVIFYENK                              5029500      4781000      5784600
      IQQGTDLAEVAPILCAGVTVYK                8426200      8015900      7577400
      IVGLSELPK                            18468000     21092000     20688000
      NMVSDIQEATK                           5963000      6346800      4587900
      VLGIDAGEEK                                 NA           NA           NA
                                       00050amol_R1 00050amol_R2 00050amol_R3
      AAADALSDLEIK                         15959000           NA           NA
      AAADALSDLEIKDSK                       1692400           NA           NA
      AEWALR                                4391400      5061500      4888800
      DEGLHTDFACLLFAHLK                    10591000      7086500      5466000
      DIHDWNNR                                   NA      1863500      1890800
      ELETLREENR                            6494400      7013800      6783800
      ESEFLFNAIHTIPEIGEK                   29300000     30287000     25566000
      GMMPGLTFSNELICR                       4700600           NA           NA
      IVTEAVEIEQR                          23396000     21146000     19843000
      LLVAFGNK                              6210100      5933900      6785700
      LLVAFGNKK                             1805500      2323300      3048700
      NKPDPAIVEK                           17325000     15126000     19103000
      TNFFEK                                4757300      6086000      5267300
      TVLFPIK                               8121400      8950100     10062000
      VENPFDFMENISLAGK                     11243000     12737000     13617000
      WIQDADALFGER                         13179000     12271000     13676000
      YFLDALPVALLGMNADLMNQYVEFVADR         25257000     26369000     25570000
      AANLGGVAVSGLEMAQNSQK                       NA      9790800      7712000
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                           24702000     21271000     24710000
      FHPSVNLSILK                          10793000      7127700      7376100
      FLGFEQIFK                            50225000     51168000     50125000
      GANIASFVMVADAMLDQGDVF               101090000     82319000     79792000
      GCIISETGITSEQIHDIASAK                 8236100      9479200      7203200
      GGLCVDLK                              9623000      9412400     10337000
      ICYAFMR                               8463300      9066000      9133000
      NSWEGVLTGK                            9933800     12459000     13323000
      SLEEIVDEYSTFSESK                      8531500      7086900      5533000
      VLPIVSVPER                           27842000     28190000     31250000
      VTISGSGNVAQYAALK                      3415800      2767600           NA
      VTWENDNGEQEVAQGYR                     2741900      2620000      2972900
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK      7305500      5124500      4483200
      AANLGGVAVSGLEMAQNSQR                 24970000     25107000     27161000
      ALVAQGVK                                   NA       627040           NA
      FIAEGSNMGSTPEAIAVFETAR               33815000     36216000     36718000
      GANIASFIK                            13626000     15858000     16235000
      GCIISETGITSEQVADISSAK                 6216600      4625600           NA
      HIGQDTDVPAGDIGVGGR                   36853000     33019000     31281000
      IMINCFNECIDYAK                       10727000     12729000     12141000
      ITWTSER                                    NA     16764000     15870000
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        77437000           NA     91002000
      SLEQIVNEYSTFSENK                     41118000      8652100     41947000
      STATGPSEAVWYGPPK                     35420000     34634000     35493000
      VDIALPCATQNEVSGEEAK                   6883200     40397000     42828000
      VIELGGTVVSLSDSK                      15966000     13700000     12572000
      VQYIAGARPWTHVQK                            NA           NA           NA
      VTWENDKGEQEVAQGYR                    49706000     45450000     45511000
      AAGLTAAYAR                           42384000     51513000     50661000
      APEAEQVLSAAATFPIAQPATDVEAR            4065200      2493700      2703800
      AVQDNGESAFR                          13853000     15309000     13296000
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        25646000     20835000     27629000
      GFTLAEVK                             26552000     31188000     31481000
      IAPRPLDLLRPVVR                       16224000      9072200      9069600
      IIVFPR                                1713200      1827500     39743000
      NQEIFDANVQR                          62089000     71286000     72922000
      TIGIAVDHR                            27022000     27046000     23572000
      VHFDQAGK                             10752000     13184000     13631000
      VHFDQAGKK                             2167800      1738600      1553600
      ANELLINVK                           141310000    148290000    156950000
      ANGTTVLVGMPAGAK                      17328000     17065000     17900000
      ATDGGAHGVINVSVSEAAIEASTR             95165000     83763000     77765000
      CCSDVFNQVVK                          54871000     55217000     59297000
      DIVGAVLK                             45878000     50659000     50399000
      EALDFFAR                            173510000    174430000    182280000
      EKDIVGAVLK                           47700000     29396000     29120000
      GVIFYESHGK                          175120000    188140000    146880000
      IGDYAGIK                                   NA    107540000    110400000
      LPLVGGHEGAGVVVGMGENVK                97824000     91234000    117110000
      SANLMAGHWVAISGAAGGLGSLAVQYAK         47157000     38196000     38193000
      SIGGEVFIDFTK                         68462000     36787000     55916000
      SIPETQK                              25681000     26642000     51120000
      SISIVGSYVGNR                        228880000    220110000    233240000
      VLGIDGGEGK                            8615500      7892100      8450300
      VLGIDGGEGKEELFR                     346850000    381800000    418900000
      VVGLSTLPEIYEK                       430530000    436840000    478320000
      YSGVCHTDLHAWHGDWPLPVK                64912000    102920000     91098000
      ANGTVVLVGLPAGAK                       4806300      4870300      4820400
      CSSDVFNHVVK                           5009000      3682300      4707200
      DIPVPKPKPNELLINVK                     6549800      5467100      3854200
      VVGLSSLPEIYEK                        26152000     26500000     30764000
      DIPVPEPKPNEILINVK                     3252300     23229000     18069000
      EALDFFSR                              5981300      5778900      5265400
      GVIFYENK                              5578200      6298300      5816500
      IQQGTDLAEVAPILCAGVTVYK                5917100      7120900      7280200
      IVGLSELPK                            12722000     14532000     14635000
      NMVSDIQEATK                           5816400      6432300      6460400
      VLGIDAGEEK                                 NA           NA           NA

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
                                                                                             Proteins
      AAADALSDLEIK                                                               sp|P09938|RIR2_YEAST
      AAADALSDLEIKDSK                                                            sp|P09938|RIR2_YEAST
      AEWALR                                                                     sp|P09938|RIR2_YEAST
      DEGLHTDFACLLFAHLK                                                          sp|P09938|RIR2_YEAST
      DIHDWNNR                                                                   sp|P09938|RIR2_YEAST
      ELETLREENR                                                                 sp|P09938|RIR2_YEAST
      ESEFLFNAIHTIPEIGEK                                                         sp|P09938|RIR2_YEAST
      GMMPGLTFSNELICR                                                            sp|P09938|RIR2_YEAST
      IVTEAVEIEQR                                                                sp|P09938|RIR2_YEAST
      LLVAFGNK                                                                   sp|P09938|RIR2_YEAST
      LLVAFGNKK                                                                  sp|P09938|RIR2_YEAST
      NKPDPAIVEK                                                                 sp|P09938|RIR2_YEAST
      TNFFEK                                                                     sp|P09938|RIR2_YEAST
      TVLFPIK                                                                    sp|P09938|RIR2_YEAST
      VENPFDFMENISLAGK                                                           sp|P09938|RIR2_YEAST
      WIQDADALFGER                                                               sp|P09938|RIR2_YEAST
      YFLDALPVALLGMNADLMNQYVEFVADR                                               sp|P09938|RIR2_YEAST
      AANLGGVAVSGLEMAQNSQK                                                       sp|P39708|DHE5_YEAST
      DAVWFGPPK                                                                  sp|P39708|DHE5_YEAST
      EIGYLFGAYR                                            sp|P07262|DHE4_YEAST;sp|P39708|DHE5_YEAST
      FHPSVNLSILK                                           sp|P07262|DHE4_YEAST;sp|P39708|DHE5_YEAST
      FLGFEQIFK                                             sp|P07262|DHE4_YEAST;sp|P39708|DHE5_YEAST
      GANIASFVMVADAMLDQGDVF                                                      sp|P39708|DHE5_YEAST
      GCIISETGITSEQIHDIASAK                                                      sp|P39708|DHE5_YEAST
      GGLCVDLK                                              sp|P07262|DHE4_YEAST;sp|P39708|DHE5_YEAST
      ICYAFMR                                               sp|P07262|DHE4_YEAST;sp|P39708|DHE5_YEAST
      NSWEGVLTGK                                            sp|P07262|DHE4_YEAST;sp|P39708|DHE5_YEAST
      SLEEIVDEYSTFSESK                                                           sp|P39708|DHE5_YEAST
      VLPIVSVPER                                            sp|P07262|DHE4_YEAST;sp|P39708|DHE5_YEAST
      VTISGSGNVAQYAALK                                      sp|P07262|DHE4_YEAST;sp|P39708|DHE5_YEAST
      VTWENDNGEQEVAQGYR                                                          sp|P39708|DHE5_YEAST
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK                                           sp|P39708|DHE5_YEAST
      AANLGGVAVSGLEMAQNSQR                                                       sp|P07262|DHE4_YEAST
      ALVAQGVK                                                                   sp|P07262|DHE4_YEAST
      FIAEGSNMGSTPEAIAVFETAR                                                     sp|P07262|DHE4_YEAST
      GANIASFIK                                                                  sp|P07262|DHE4_YEAST
      GCIISETGITSEQVADISSAK                                                      sp|P07262|DHE4_YEAST
      HIGQDTDVPAGDIGVGGR                                                         sp|P07262|DHE4_YEAST
      IMINCFNECIDYAK                                                             sp|P07262|DHE4_YEAST
      ITWTSER                                                                    sp|P07262|DHE4_YEAST
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR                                              sp|P07262|DHE4_YEAST
      SLEQIVNEYSTFSENK                                                           sp|P07262|DHE4_YEAST
      STATGPSEAVWYGPPK                                                           sp|P07262|DHE4_YEAST
      VDIALPCATQNEVSGEEAK                                                        sp|P07262|DHE4_YEAST
      VIELGGTVVSLSDSK                                                            sp|P07262|DHE4_YEAST
      VQYIAGARPWTHVQK                                                            sp|P07262|DHE4_YEAST
      VTWENDKGEQEVAQGYR                                                          sp|P07262|DHE4_YEAST
      AAGLTAAYAR                                          sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST
      APEAEQVLSAAATFPIAQPATDVEAR                          sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST
      AVQDNGESAFR                                         sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST
      DGKAPEAEQVLSAAATFPIAQPATDVEAR                                             sp|P40212|RL13B_YEAST
      GFTLAEVK                                            sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST
      IAPRPLDLLRPVVR                                      sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST
      IIVFPR                                              sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST
      NQEIFDANVQR                                         sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST
      TIGIAVDHR                                           sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST
      VHFDQAGK                                            sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST
      VHFDQAGKK                                           sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST
      ANELLINVK                                                                  sp|P00330|ADH1_YEAST
      ANGTTVLVGMPAGAK                                                            sp|P00330|ADH1_YEAST
      ATDGGAHGVINVSVSEAAIEASTR                                                   sp|P00330|ADH1_YEAST
      CCSDVFNQVVK                                                                sp|P00330|ADH1_YEAST
      DIVGAVLK                                                                   sp|P00330|ADH1_YEAST
      EALDFFAR                         sp|P00330|ADH1_YEAST;sp|P38113|ADH5_YEAST;sp|P00331|ADH2_YEAST
      EKDIVGAVLK                                                                 sp|P00330|ADH1_YEAST
      GVIFYESHGK                                                                 sp|P00330|ADH1_YEAST
      IGDYAGIK                                              sp|P00330|ADH1_YEAST;sp|P00331|ADH2_YEAST
      LPLVGGHEGAGVVVGMGENVK                                 sp|P00330|ADH1_YEAST;sp|P00331|ADH2_YEAST
      SANLMAGHWVAISGAAGGLGSLAVQYAK                                               sp|P00330|ADH1_YEAST
      SIGGEVFIDFTK                                                               sp|P00330|ADH1_YEAST
      SIPETQK                                               sp|P00330|ADH1_YEAST;sp|P00331|ADH2_YEAST
      SISIVGSYVGNR                                          sp|P00330|ADH1_YEAST;sp|P00331|ADH2_YEAST
      VLGIDGGEGK                                                                 sp|P00330|ADH1_YEAST
      VLGIDGGEGKEELFR                                                            sp|P00330|ADH1_YEAST
      VVGLSTLPEIYEK                                                              sp|P00330|ADH1_YEAST
      YSGVCHTDLHAWHGDWPLPVK                                 sp|P00330|ADH1_YEAST;sp|P07246|ADH3_YEAST
      ANGTVVLVGLPAGAK                                                            sp|P00331|ADH2_YEAST
      CSSDVFNHVVK                                                                sp|P00331|ADH2_YEAST
      DIPVPKPKPNELLINVK                                                          sp|P00331|ADH2_YEAST
      VVGLSSLPEIYEK                                                              sp|P00331|ADH2_YEAST
      DIPVPEPKPNEILINVK                                                          sp|P07246|ADH3_YEAST
      EALDFFSR                                                                   sp|P07246|ADH3_YEAST
      GVIFYENK                                                                   sp|P07246|ADH3_YEAST
      IQQGTDLAEVAPILCAGVTVYK                                                     sp|P07246|ADH3_YEAST
      IVGLSELPK                                                                  sp|P07246|ADH3_YEAST
      NMVSDIQEATK                                                                sp|P07246|ADH3_YEAST
      VLGIDAGEEK                                                                 sp|P07246|ADH3_YEAST
                                         Score
      AAADALSDLEIK                     201.750
      AAADALSDLEIKDSK                   84.297
      AEWALR                           109.700
      DEGLHTDFACLLFAHLK                163.620
      DIHDWNNR                          47.574
      ELETLREENR                       167.230
      ESEFLFNAIHTIPEIGEK               232.130
      GMMPGLTFSNELICR                  115.570
      IVTEAVEIEQR                      188.910
      LLVAFGNK                         100.070
      LLVAFGNKK                         78.616
      NKPDPAIVEK                        71.451
      TNFFEK                           124.340
      TVLFPIK                           92.191
      VENPFDFMENISLAGK                 217.190
      WIQDADALFGER                     199.800
      YFLDALPVALLGMNADLMNQYVEFVADR     110.310
      AANLGGVAVSGLEMAQNSQK             211.490
      DAVWFGPPK                         73.665
      EIGYLFGAYR                       161.020
      FHPSVNLSILK                      208.330
      FLGFEQIFK                        119.210
      GANIASFVMVADAMLDQGDVF             32.565
      GCIISETGITSEQIHDIASAK             94.544
      GGLCVDLK                         122.130
      ICYAFMR                          126.050
      NSWEGVLTGK                       155.110
      SLEEIVDEYSTFSESK                 284.150
      VLPIVSVPER                       120.150
      VTISGSGNVAQYAALK                 171.620
      VTWENDNGEQEVAQGYR                174.940
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK  27.863
      AANLGGVAVSGLEMAQNSQR             342.280
      ALVAQGVK                         154.090
      FIAEGSNMGSTPEAIAVFETAR           309.150
      GANIASFIK                        113.180
      GCIISETGITSEQVADISSAK            328.870
      HIGQDTDVPAGDIGVGGR               247.170
      IMINCFNECIDYAK                   220.290
      ITWTSER                          107.590
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR    104.260
      SLEQIVNEYSTFSENK                 383.440
      STATGPSEAVWYGPPK                 190.360
      VDIALPCATQNEVSGEEAK              242.010
      VIELGGTVVSLSDSK                  273.150
      VQYIAGARPWTHVQK                   51.211
      VTWENDKGEQEVAQGYR                256.130
      AAGLTAAYAR                       158.070
      APEAEQVLSAAATFPIAQPATDVEAR        83.221
      AVQDNGESAFR                      154.840
      DGKAPEAEQVLSAAATFPIAQPATDVEAR    115.790
      GFTLAEVK                         122.740
      IAPRPLDLLRPVVR                    90.793
      IIVFPR                            98.352
      NQEIFDANVQR                      225.010
      TIGIAVDHR                        141.520
      VHFDQAGK                          84.916
      VHFDQAGKK                        105.100
      ANELLINVK                        203.350
      ANGTTVLVGMPAGAK                  176.870
      ATDGGAHGVINVSVSEAAIEASTR         200.480
      CCSDVFNQVVK                      282.510
      DIVGAVLK                         124.340
      EALDFFAR                         165.110
      EKDIVGAVLK                       156.350
      GVIFYESHGK                       164.960
      IGDYAGIK                         114.780
      LPLVGGHEGAGVVVGMGENVK            139.120
      SANLMAGHWVAISGAAGGLGSLAVQYAK      89.979
      SIGGEVFIDFTK                     150.460
      SIPETQK                           72.138
      SISIVGSYVGNR                     280.960
      VLGIDGGEGK                        80.455
      VLGIDGGEGKEELFR                  193.150
      VVGLSTLPEIYEK                    203.020
      YSGVCHTDLHAWHGDWPLPVK            200.900
      ANGTVVLVGLPAGAK                   94.538
      CSSDVFNHVVK                       98.040
      DIPVPKPKPNELLINVK                 52.372
      VVGLSSLPEIYEK                    202.560
      DIPVPEPKPNEILINVK                124.380
      EALDFFSR                          85.064
      GVIFYENK                          96.190
      IQQGTDLAEVAPILCAGVTVYK           226.740
      IVGLSELPK                        110.120
      NMVSDIQEATK                       89.679
      VLGIDAGEEK                        59.350

---

    Code
      as.data.frame(tail(SummarizedExperiment::colData(D2), n = 1000))
    Output
                         sample
      12500amol_R1 12500amol_R1
      12500amol_R2 12500amol_R2
      12500amol_R3 12500amol_R3
      00125amol_R1 00125amol_R1
      00125amol_R2 00125amol_R2
      00125amol_R3 00125amol_R3
      25000amol_R1 25000amol_R1
      25000amol_R2 25000amol_R2
      25000amol_R3 25000amol_R3
      02500amol_R1 02500amol_R1
      02500amol_R2 02500amol_R2
      02500amol_R3 02500amol_R3
      00250amol_R1 00250amol_R1
      00250amol_R2 00250amol_R2
      00250amol_R3 00250amol_R3
      50000amol_R1 50000amol_R1
      50000amol_R2 50000amol_R2
      50000amol_R3 50000amol_R3
      05000amol_R1 05000amol_R1
      05000amol_R2 05000amol_R2
      05000amol_R3 05000amol_R3
      00500amol_R1 00500amol_R1
      00500amol_R2 00500amol_R2
      00500amol_R3 00500amol_R3
      00050amol_R1 00050amol_R1
      00050amol_R2 00050amol_R2
      00050amol_R3 00050amol_R3

