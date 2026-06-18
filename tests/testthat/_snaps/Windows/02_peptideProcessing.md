# normalize peptide data

    Code
      SummarizedExperiment::assays(D_norm_loess)$intensities_norm
    Output
                                       12500amol_R1 12500amol_R2 12500amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                       7489325      7607065    5367203.1
      AEWALR                                4767329      4686172    5343109.4
      DEGLHTDFACLLFAHLK                     4215870      3593719     533167.4
      DIHDWNNR                              1725574      1845995           NA
      ELETLREENR                                 NA           NA           NA
      ESEFLFNAIHTIPEIGEK                   32444583     34689550   33017820.0
      GMMPGLTFSNELICR                       8548691      7109225    8394714.9
      IVTEAVEIEQR                          18145994     17973084   17855343.0
      LLVAFGNK                             11201984     12130294   11954322.4
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                           16026462     15298661   16036457.1
      TNFFEK                                6643283      7133204    5256534.5
      TVLFPIK                              13309919     12642767   12478198.1
      VENPFDFMENISLAGK                     16293262     16989631   17986154.4
      WIQDADALFGER                         17172467     17626754   17009773.9
      YFLDALPVALLGMNADLMNQYVEFVADR         25379727     28588104   27282033.5
      AANLGGVAVSGLEMAQNSQK                  8940105      3404003   10139854.2
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           4136386           NA    5043326.6
      FLGFEQIFK                            39822973     39787388   41160417.3
      GANIASFVMVADAMLDQGDVF                56967887     56562859   57990635.2
      GCIISETGITSEQIHDIASAK                 8072626      8067019    8616572.8
      GGLCVDLK                                   NA      9705722           NA
      ICYAFMR                               8990873      9041759    9738469.2
      NSWEGVLTGK                           11681321     11140904   13598402.6
      SLEEIVDEYSTFSESK                      8341678      8189566   10497415.7
      VLPIVSVPER                           39106436     36122874   38924830.3
      VTISGSGNVAQYAALK                      4295351       138558    6029471.2
      VTWENDNGEQEVAQGYR                     3268490           NA    3614412.1
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 28231351     30229200   26794481.8
      ALVAQGVK                             10591041     12363150   12020035.3
      FIAEGSNMGSTPEAIAVFETAR               22809430     39806035   43540097.1
      GANIASFIK                            18607205     14629885   15766046.4
      GCIISETGITSEQVADISSAK                 3069545      3721888    4188169.2
      HIGQDTDVPAGDIGVGGR                   27113789     31182405   33465304.9
      IMINCFNECIDYAK                       13572131     12911574   11946588.9
      ITWTSER                              19535232     16303079   15359188.9
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        10974842     89914206   94021147.0
      SLEQIVNEYSTFSENK                     53726038     69000561   57274857.5
      STATGPSEAVWYGPPK                     48375353     39360472   40016147.0
      VDIALPCATQNEVSGEEAK                  47868040     45162249   50535809.3
      VIELGGTVVSLSDSK                      13158834     13704038   15008507.7
      VQYIAGARPWTHVQK                            NA      6675490    4389794.4
      VTWENDKGEQEVAQGYR                    33174921     31883118   37891241.6
      AAGLTAAYAR                           63911309     61676209   53246397.7
      APEAEQVLSAAATFPIAQPATDVEAR           28312199     27229128   26468091.8
      AVQDNGESAFR                          10160106     12429942   12071222.7
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        26744962     26915338   26999609.3
      GFTLAEVK                             27652223     28430464   28320749.2
      IAPRPLDLLRPVVR                       16819023     15689147   17529227.8
      IIVFPR                               48990803     46588217   49369949.1
      NQEIFDANVQR                          88192410     90637352   86843645.7
      TIGIAVDHR                            23952414     24938182   22571144.6
      VHFDQAGK                              8916916      8999322   11449450.6
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           175413197    186725233  186793096.1
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR            114334846    115317500  123264945.4
      CCSDVFNQVVK                          49955864     49865748   48172822.6
      DIVGAVLK                                   NA     52264795           NA
      EALDFFAR                            203339225    201680782  196946247.8
      EKDIVGAVLK                           37362183     36727663   41074828.4
      GVIFYESHGK                           92349150     98616437  105584801.8
      IGDYAGIK                            149453700    136361224  144632852.5
      LPLVGGHEGAGVVVGMGENVK               126319037    128210887  149749559.4
      SANLMAGHWVAISGAAGGLGSLAVQYAK         27292366     28439304   23062613.9
      SIGGEVFIDFTK                         71586865     67323643   58003469.6
      SIPETQK                              28995661     24978857   25186769.8
      SISIVGSYVGNR                        283558402    302422757  276859910.5
      VLGIDGGEGK                           20718057     19114983           NA
      VLGIDGGEGKEELFR                     320501214    330491188  337477480.3
      VVGLSTLPEIYEK                       600541629    560815056  563267041.9
      YSGVCHTDLHAWHGDWPLPVK               103325789     96693491  107190860.3
      ANGTVVLVGLPAGAK                       4761599      4460463    4664558.8
      CSSDVFNHVVK                                NA           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        38321275     35122040   37414397.1
      DIPVPEPKPNEILINVK                    17911208     19531958   15745072.5
      EALDFFSR                              7074888      5794466    6418783.6
      GVIFYENK                              5169306      5250695    5335378.9
      IQQGTDLAEVAPILCAGVTVYK                7802892      7887646    6024976.5
      IVGLSELPK                            16283393     15286875   17265093.8
      NMVSDIQEATK                           5406573      6758562    5853129.8
      VLGIDAGEEK                                 NA           NA           NA
                                       125amol_R1  125amol_R2  125amol_R3
      AAADALSDLEIK                             NA          NA          NA
      AAADALSDLEIKDSK                    13041901  12978564.8  13309129.3
      AEWALR                              6053901   5616839.2   5224349.6
      DEGLHTDFACLLFAHLK                   2970527   4584935.8   4029347.1
      DIHDWNNR                            2000377          NA   2340378.7
      ELETLREENR                          5154088   4621965.7   6482132.2
      ESEFLFNAIHTIPEIGEK                 32565395  30005324.9  37917100.1
      GMMPGLTFSNELICR                     6574659   7337990.7   7858718.4
      IVTEAVEIEQR                        16812280  17927200.5  17833816.2
      LLVAFGNK                            6363777   7211833.9   7519148.5
      LLVAFGNKK                           1597502   1782824.1   1663657.3
      NKPDPAIVEK                         18429313  17556529.7  16319503.3
      TNFFEK                              5832809   4987826.8   4938628.5
      TVLFPIK                            11092844  10105316.6  10142811.3
      VENPFDFMENISLAGK                   17764469  16456256.3  14602843.7
      WIQDADALFGER                       15830177  15560010.1  16277209.6
      YFLDALPVALLGMNADLMNQYVEFVADR       24430594  21059770.1  20295708.4
      AANLGGVAVSGLEMAQNSQK                9785062   9176426.4   9593603.7
      DAVWFGPPK                                NA          NA          NA
      EIGYLFGAYR                         23815029  24103719.1  22765845.7
      FHPSVNLSILK                         4933063   4664117.0   5066383.9
      FLGFEQIFK                          49255542  43784814.0  45082635.4
      GANIASFVMVADAMLDQGDVF              83259687  83726928.9  76969370.3
      GCIISETGITSEQIHDIASAK               8335613   7374872.1   8590504.7
      GGLCVDLK                           10565237          NA          NA
      ICYAFMR                            11427237  12962643.3  11619700.5
      NSWEGVLTGK                         12520442  17072173.8          NA
      SLEEIVDEYSTFSESK                    9278312   9553986.6   8930342.9
      VLPIVSVPER                         32868732  36800063.3  30972103.6
      VTISGSGNVAQYAALK                    6542698     99915.1          NA
      VTWENDNGEQEVAQGYR                        NA   2743143.3   2241969.3
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK         NA          NA   1784617.2
      AANLGGVAVSGLEMAQNSQR                9330629  33000616.5  32744214.8
      ALVAQGVK                           13850234  13631875.0  13432524.3
      FIAEGSNMGSTPEAIAVFETAR             40036784  39940956.3  41437759.0
      GANIASFIK                          17833728  17770883.4  18599309.2
      GCIISETGITSEQVADISSAK               2894457   2901773.9   2636308.2
      HIGQDTDVPAGDIGVGGR                 33669234  33422622.6  31768071.5
      IMINCFNECIDYAK                     12009153  12457684.6  12002011.2
      ITWTSER                            16555098  16239949.1    508059.5
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR       9885907 112578972.5  96390092.7
      SLEQIVNEYSTFSENK                   46993164  50182398.5  43233404.2
      STATGPSEAVWYGPPK                   38141064  41517999.2  44792845.2
      VDIALPCATQNEVSGEEAK                47403753  50829878.8  47757261.1
      VIELGGTVVSLSDSK                    12538825  12635217.2  12660385.5
      VQYIAGARPWTHVQK                     3735310          NA   3754307.0
      VTWENDKGEQEVAQGYR                  42399742  39115675.5  42568405.7
      AAGLTAAYAR                         60750205  60057538.2  59646246.1
      APEAEQVLSAAATFPIAQPATDVEAR         18985688  23109268.2  20254406.0
      AVQDNGESAFR                        13906061  13126820.3  14213919.1
      DGKAPEAEQVLSAAATFPIAQPATDVEAR      39277840  40360805.9  39545726.3
      GFTLAEVK                           31017605  31688884.3  30555223.6
      IAPRPLDLLRPVVR                     22082765  26752474.0  26457531.6
      IIVFPR                             40002487  41568631.3   3093299.3
      NQEIFDANVQR                        83467953  80246798.4  78253819.4
      TIGIAVDHR                          19294003  18273498.3  16730563.3
      VHFDQAGK                                 NA          NA  13292733.5
      VHFDQAGKK                                NA          NA          NA
      ANELLINVK                         169723284 175577212.9 169671869.3
      ANGTTVLVGMPAGAK                          NA          NA          NA
      ATDGGAHGVINVSVSEAAIEASTR           85020998  63241687.2  79641984.1
      CCSDVFNQVVK                        57489246  58195650.1  52150932.2
      DIVGAVLK                                 NA  56615350.9  49106565.5
      EALDFFAR                          204810564 203494181.0 203031572.5
      EKDIVGAVLK                         30905412  31745998.0  33025993.1
      GVIFYESHGK                        120599767 116284171.8 107357230.5
      IGDYAGIK                          139018421 133414295.9 138430336.3
      LPLVGGHEGAGVVVGMGENVK             115912243  83723461.9 147810455.0
      SANLMAGHWVAISGAAGGLGSLAVQYAK       21882261  27264928.5  31237557.0
      SIGGEVFIDFTK                       82167898  61798965.1  72932513.3
      SIPETQK                            26170701  43853021.7  45517845.7
      SISIVGSYVGNR                      262956733 280249932.1 282635669.4
      VLGIDGGEGK                          9379426  11315222.5  10345183.8
      VLGIDGGEGKEELFR                   423289369 389844139.4 339776361.1
      VVGLSTLPEIYEK                     556942235 551383377.0 573087807.8
      YSGVCHTDLHAWHGDWPLPVK              93846753 104016482.3 115201200.3
      ANGTVVLVGLPAGAK                     4654159   5727859.3   4909999.2
      CSSDVFNHVVK                         4020076   4690232.1   4044162.6
      DIPVPKPKPNELLINVK                        NA   1171711.2          NA
      VVGLSSLPEIYEK                      35063423  34349404.5  33891256.4
      DIPVPEPKPNEILINVK                  17142402  20490423.1  19710857.9
      EALDFFSR                            6828997   7120763.1   6429381.6
      GVIFYENK                                 NA   5828471.5   5616936.9
      IQQGTDLAEVAPILCAGVTVYK              7862077   7377289.8   8061740.5
      IVGLSELPK                          17781845  18380137.7  20096261.4
      NMVSDIQEATK                         3666308   4214223.0          NA
      VLGIDAGEEK                               NA          NA          NA
                                       25000amol_R1 25000amol_R2 25000amol_R3
      AAADALSDLEIK                               NA           NA     29399758
      AAADALSDLEIKDSK                       7626932    9224596.8      5263071
      AEWALR                                4885375    5487028.2      5643348
      DEGLHTDFACLLFAHLK                     5128512    3490350.7      6255297
      DIHDWNNR                              2513906     862312.4      2635226
      ELETLREENR                                 NA    1102907.5           NA
      ESEFLFNAIHTIPEIGEK                   32971637   27401873.2     34415540
      GMMPGLTFSNELICR                       7991037    8040233.8      6988032
      IVTEAVEIEQR                          18289328   19694045.3     20950238
      LLVAFGNK                             10741688   13603615.8     11137428
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                           14675556   15058780.5     17922495
      TNFFEK                                5841404    6461503.0      7456515
      TVLFPIK                              11839020   13157318.7     12344145
      VENPFDFMENISLAGK                     17053571   15159649.9     13535800
      WIQDADALFGER                         18009980   16757166.1     15732065
      YFLDALPVALLGMNADLMNQYVEFVADR         29990818   21028134.7     21416232
      AANLGGVAVSGLEMAQNSQK                  9066051    9529002.4      7960614
      DAVWFGPPK                           178406359           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                           6012713    7066933.9      6685306
      FLGFEQIFK                            38662626   38069595.9     34190506
      GANIASFVMVADAMLDQGDVF                55934060   43856182.8     44954780
      GCIISETGITSEQIHDIASAK                 7723456    6362956.7      8206230
      GGLCVDLK                              9913917   10461559.5     10960814
      ICYAFMR                               8449764    6076897.6      8379030
      NSWEGVLTGK                           13368976   13600805.8     12164402
      SLEEIVDEYSTFSESK                     10037302    6699082.3      8321248
      VLPIVSVPER                           40699541   49701568.8     40619431
      VTISGSGNVAQYAALK                      1836578    7264875.5      1299357
      VTWENDNGEQEVAQGYR                     2504429           NA      2246853
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 30613549   29632186.6     26490583
      ALVAQGVK                             12051660   15689995.6     15884181
      FIAEGSNMGSTPEAIAVFETAR               39509947   41586767.7     36745978
      GANIASFIK                            16928338   16087826.3           NA
      GCIISETGITSEQVADISSAK                 3029898    3440416.9      3739464
      HIGQDTDVPAGDIGVGGR                   30090328   32667481.2     35058807
      IMINCFNECIDYAK                       11024925   10729186.3     10339968
      ITWTSER                              16775241   20009946.4     18719058
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR       102833320   71456033.2     84776072
      SLEQIVNEYSTFSENK                     57524621   51591898.8     56664936
      STATGPSEAVWYGPPK                     43424490   50440366.2     41087636
      VDIALPCATQNEVSGEEAK                  51935120   54422729.7     49593319
      VIELGGTVVSLSDSK                      13284660   14522261.5     13446024
      VQYIAGARPWTHVQK                       3670269    6441393.5      2885840
      VTWENDKGEQEVAQGYR                    34748566   29803621.5     31248285
      AAGLTAAYAR                           60161935   69116949.9     68838349
      APEAEQVLSAAATFPIAQPATDVEAR           30030125   27609400.5     29549648
      AVQDNGESAFR                          11371175   12422974.2     13324207
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        18655039   22111884.7     22209354
      GFTLAEVK                             29684530   27065974.1     26950013
      IAPRPLDLLRPVVR                       16298683   13039857.2     16988127
      IIVFPR                               51967257   57576567.5     51267475
      NQEIFDANVQR                          84293428   92008299.3     89682178
      TIGIAVDHR                            24269675   23698744.9     26543444
      VHFDQAGK                             11458054    9909970.1     11178755
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           190021008  198794074.4    187696553
      ANGTTVLVGMPAGAK                            NA           NA           NA
      ATDGGAHGVINVSVSEAAIEASTR             86820670   91087906.7    114548484
      CCSDVFNQVVK                          48219067   52005351.9     48157129
      DIVGAVLK                             50424573   58623866.8     50947751
      EALDFFAR                            203717921  225364465.5    190002341
      EKDIVGAVLK                           38712168   30732198.2     31229165
      GVIFYESHGK                           99145280   96028911.5     99916594
      IGDYAGIK                            152026623  180641399.6    177008247
      LPLVGGHEGAGVVVGMGENVK               130117863  100044341.5    120808397
      SANLMAGHWVAISGAAGGLGSLAVQYAK         27711555   41449554.9     31301137
      SIGGEVFIDFTK                         95608644   58917745.4     50296135
      SIPETQK                              23239585   19862545.0     13343614
      SISIVGSYVGNR                        293449742  303964576.0    297094787
      VLGIDGGEGK                           19362185   26765985.6     26679050
      VLGIDGGEGKEELFR                     327495458  287038370.1    299346864
      VVGLSTLPEIYEK                       544538633  585909996.6    562733015
      YSGVCHTDLHAWHGDWPLPVK               117704298   91041000.3    110941834
      ANGTVVLVGLPAGAK                       4458327    5254372.4      4961309
      CSSDVFNHVVK                           4290928           NA           NA
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        40409941   40702790.7     38110287
      DIPVPEPKPNEILINVK                    20495736   21901324.8     20333057
      EALDFFSR                              6476848    6722983.5      5024664
      GVIFYENK                              4703566    4313196.9      6231937
      IQQGTDLAEVAPILCAGVTVYK                8188956    6823295.3      5978849
      IVGLSELPK                            18047578   16452086.8     15075144
      NMVSDIQEATK                           5871209    6627152.0      6292980
      VLGIDAGEEK                                 NA           NA           NA
                                       2500amol_R1 2500amol_R2 2500amol_R3
      AAADALSDLEIK                              NA          NA          NA
      AAADALSDLEIKDSK                    5714895.0   7412859.9     8323247
      AEWALR                             5437371.2   4759844.2     4687328
      DEGLHTDFACLLFAHLK                  4704902.3   5284551.7     3701090
      DIHDWNNR                            528414.8          NA          NA
      ELETLREENR                                NA          NA     4475112
      ESEFLFNAIHTIPEIGEK                33332941.5  26526564.1    35233678
      GMMPGLTFSNELICR                    8372378.6   9351807.6     8949268
      IVTEAVEIEQR                       16788306.6  18097768.5    18150568
      LLVAFGNK                          10772601.0  10064115.8     9102412
      LLVAFGNKK                           640883.2          NA          NA
      NKPDPAIVEK                        13862242.7  14399683.1    17719191
      TNFFEK                             6550470.0   5755964.4     5927178
      TVLFPIK                           10153058.2  10821968.1    11671158
      VENPFDFMENISLAGK                  16247278.9  15259613.4    15813557
      WIQDADALFGER                      18552752.9  21121961.6    18732935
      YFLDALPVALLGMNADLMNQYVEFVADR      23560114.2  23763038.6    24349589
      AANLGGVAVSGLEMAQNSQK              10753748.7  10083377.2     8032028
      DAVWFGPPK                                 NA          NA          NA
      EIGYLFGAYR                                NA          NA          NA
      FHPSVNLSILK                        5348654.6   4767786.3     7520494
      FLGFEQIFK                         40325799.5  43302314.7    39014591
      GANIASFVMVADAMLDQGDVF             53607086.1  58725440.3    56952391
      GCIISETGITSEQIHDIASAK              6091271.1   7855327.2     8234518
      GGLCVDLK                           9745216.7   8917139.6    10735804
      ICYAFMR                            8460053.0   8407758.2     9297939
      NSWEGVLTGK                        12557775.1  12602459.6    11338263
      SLEEIVDEYSTFSESK                   8391215.9  10625759.0     7975351
      VLPIVSVPER                        40861208.7  40311373.1    42622995
      VTISGSGNVAQYAALK                   5698308.8    208949.8          NA
      VTWENDNGEQEVAQGYR                  3067734.5   3169706.0     3137459
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK          NA          NA          NA
      AANLGGVAVSGLEMAQNSQR              30482840.8  34398601.1    32481849
      ALVAQGVK                                  NA  11665060.4    12797820
      FIAEGSNMGSTPEAIAVFETAR            35695567.8  38067089.6    35966064
      GANIASFIK                         17189545.0  17391342.6    18097550
      GCIISETGITSEQVADISSAK              3489029.9   3402323.8     3508332
      HIGQDTDVPAGDIGVGGR                33190618.6  28895299.5    31876379
      IMINCFNECIDYAK                    12842923.8  12913853.4    13063854
      ITWTSER                                   NA  16237047.4    16277578
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     94357104.7 100066631.6    84033147
      SLEQIVNEYSTFSENK                  54342353.3  54131516.5    58382151
      STATGPSEAVWYGPPK                  49572961.3  50831826.3    44554442
      VDIALPCATQNEVSGEEAK               54775630.2  58008598.0    53904785
      VIELGGTVVSLSDSK                   13569105.6  13383549.4    13192828
      VQYIAGARPWTHVQK                    8043263.2          NA     6126999
      VTWENDKGEQEVAQGYR                 33702820.3  32741866.8    33647535
      AAGLTAAYAR                        57871818.3  57524482.3    55836525
      APEAEQVLSAAATFPIAQPATDVEAR        25156320.4  27827693.3    28493609
      AVQDNGESAFR                       13176609.9  12197062.7    12968943
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     26725376.7  22818986.1    20861009
      GFTLAEVK                          28712916.7  32088546.4    29762140
      IAPRPLDLLRPVVR                     8434923.0  19215713.6     9966184
      IIVFPR                            57941020.3  57304266.8    51163569
      NQEIFDANVQR                       93293911.4  91871742.3    88074549
      TIGIAVDHR                         21072371.5  23117220.1    22998159
      VHFDQAGK                           9085825.9          NA    10352931
      VHFDQAGKK                                 NA          NA          NA
      ANELLINVK                        171750516.0 157142993.7   160394454
      ANGTTVLVGMPAGAK                           NA          NA          NA
      ATDGGAHGVINVSVSEAAIEASTR          77963789.4  97835410.8    91155116
      CCSDVFNQVVK                       54167513.0  53156865.0    52967568
      DIVGAVLK                          49305867.7          NA    48835095
      EALDFFAR                         217525315.2 213397782.4   206540360
      EKDIVGAVLK                        28878462.0  31112090.1    36104239
      GVIFYESHGK                        93466082.3  91952822.9    93547069
      IGDYAGIK                         151454908.3 137397461.0   147332546
      LPLVGGHEGAGVVVGMGENVK            111484138.4 122828404.2   157536907
      SANLMAGHWVAISGAAGGLGSLAVQYAK      35680599.8  21499696.1    28371037
      SIGGEVFIDFTK                      55180112.4  74070680.4    72810776
      SIPETQK                           37891046.4  21533145.2    23561713
      SISIVGSYVGNR                     295639177.6 286424600.5   272751958
      VLGIDGGEGK                        18817229.8  18559786.0    20612400
      VLGIDGGEGKEELFR                  315295519.4 336889630.0   356117633
      VVGLSTLPEIYEK                    603730851.5 579031244.8   553723465
      YSGVCHTDLHAWHGDWPLPVK            161170999.7 160911205.0   145735496
      ANGTVVLVGLPAGAK                    4770901.5   4045025.2     4975299
      CSSDVFNHVVK                               NA          NA          NA
      DIPVPKPKPNELLINVK                         NA          NA          NA
      VVGLSSLPEIYEK                     35681845.2  38823267.1    34429563
      DIPVPEPKPNEILINVK                 20998288.9  22970744.1    21707069
      EALDFFSR                           7059423.1   6330589.1     8580457
      GVIFYENK                           4823890.6   5231286.9     5525557
      IQQGTDLAEVAPILCAGVTVYK             6622524.2   8248172.8     6752167
      IVGLSELPK                         18271273.2  19207124.3    19380364
      NMVSDIQEATK                        6112770.6   6829729.2     5958651
      VLGIDAGEEK                                NA          NA          NA
                                        250amol_R1 250amol_R2  250amol_R3
      AAADALSDLEIK                              NA         NA          NA
      AAADALSDLEIKDSK                   10396291.2   13552580  12692279.6
      AEWALR                             4991268.5    6014079   5034461.6
      DEGLHTDFACLLFAHLK                         NA         NA          NA
      DIHDWNNR                                  NA    1786889   1870160.6
      ELETLREENR                         4006530.5    3929205   4014692.5
      ESEFLFNAIHTIPEIGEK                27357049.5   24069599  35160323.2
      GMMPGLTFSNELICR                    7458944.9         NA   7786589.8
      IVTEAVEIEQR                       18458574.0   16776054  18614311.9
      LLVAFGNK                           6876954.0    6260006   7383849.6
      LLVAFGNKK                                 NA    1698914          NA
      NKPDPAIVEK                        16794050.3   19286735  18225336.4
      TNFFEK                             6519141.7    5625847   5487528.7
      TVLFPIK                            9228187.4    9550729  11626352.6
      VENPFDFMENISLAGK                  17729458.4   17447192  17752585.9
      WIQDADALFGER                      14718856.5   14539567  16692310.0
      YFLDALPVALLGMNADLMNQYVEFVADR      24260410.0   22784190  23067157.3
      AANLGGVAVSGLEMAQNSQK              10586905.8   10157127   6058514.0
      DAVWFGPPK                                 NA         NA          NA
      EIGYLFGAYR                                NA   20427121  20668569.9
      FHPSVNLSILK                        4709423.1    6464082   5413956.1
      FLGFEQIFK                         46022191.8   45389081  43819767.5
      GANIASFVMVADAMLDQGDVF             76158914.4   86513214  80704807.0
      GCIISETGITSEQIHDIASAK              8529070.5    8741604   7134222.0
      GGLCVDLK                          11301239.0   11192254  10755319.2
      ICYAFMR                           10226042.5   10562932  12329951.5
      NSWEGVLTGK                        12164062.3   13564537  14098660.5
      SLEEIVDEYSTFSESK                   8589457.1    8113462   8339019.0
      VLPIVSVPER                        32671122.1   33118882  38575343.5
      VTISGSGNVAQYAALK                   3259856.4         NA    136915.9
      VTWENDNGEQEVAQGYR                  2859220.9    2632556   2985522.7
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK          NA         NA          NA
      AANLGGVAVSGLEMAQNSQR              36818819.3   33980012  37404323.4
      ALVAQGVK                          13673207.7   12659607  12143227.3
      FIAEGSNMGSTPEAIAVFETAR            41652398.5   27149236  41369962.5
      GANIASFIK                         18004742.2   17759614  18561254.9
      GCIISETGITSEQVADISSAK              4036246.9    3585912   3359915.6
      HIGQDTDVPAGDIGVGGR                32585427.2   36052038  30111918.8
      IMINCFNECIDYAK                    14336324.4   12950799  13098112.1
      ITWTSER                           17526846.9   16519054  16296391.7
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR    106299471.4  101646105 102361382.7
      SLEQIVNEYSTFSENK                  53657305.8   54290501  60605229.1
      STATGPSEAVWYGPPK                  45690910.6   45275429  37709868.3
      VDIALPCATQNEVSGEEAK               48129701.3   43428983  50750544.6
      VIELGGTVVSLSDSK                   13619451.3   12835112  13019261.8
      VQYIAGARPWTHVQK                   10353977.1    9945650   4594227.6
      VTWENDKGEQEVAQGYR                 35634156.4   41638101  38033386.5
      AAGLTAAYAR                        63245793.2   58963316  60990726.5
      APEAEQVLSAAATFPIAQPATDVEAR        20445100.1   17998112  14394461.6
      AVQDNGESAFR                         272322.1   13463524  11697682.9
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     27589606.1   26500602  27603120.4
      GFTLAEVK                          30328127.1   33688747  29635818.8
      IAPRPLDLLRPVVR                    13138971.4   23715390  23468948.0
      IIVFPR                            42443916.1   47390594  42003098.5
      NQEIFDANVQR                       93938256.6   85952850  82262029.8
      TIGIAVDHR                         21377812.9   22601043  19712559.6
      VHFDQAGK                          10672956.6   10594267  10342502.2
      VHFDQAGKK                                 NA         NA          NA
      ANELLINVK                        175454414.8  182719447 167837345.3
      ANGTTVLVGMPAGAK                           NA         NA          NA
      ATDGGAHGVINVSVSEAAIEASTR          81016675.8   80758141  69302619.2
      CCSDVFNQVVK                       51490927.8   57795874  56804759.5
      DIVGAVLK                          53305441.6   56211992          NA
      EALDFFAR                         205259813.3  199242927 216428051.6
      EKDIVGAVLK                        33194164.2   36146807  35738986.4
      GVIFYESHGK                       104327182.9   96151512 105268634.9
      IGDYAGIK                         142310091.3  138460951          NA
      LPLVGGHEGAGVVVGMGENVK            125297752.9  129234650  89976577.8
      SANLMAGHWVAISGAAGGLGSLAVQYAK      26648350.7   29884638  27894257.3
      SIGGEVFIDFTK                      64527497.0   62387644  55655116.8
      SIPETQK                           28231330.0   23201112  23889137.1
      SISIVGSYVGNR                     267686890.3  278248911 277746291.2
      VLGIDGGEGK                        10080194.8   10546422  10422117.8
      VLGIDGGEGKEELFR                  348786554.5  364594934 372162040.3
      VVGLSTLPEIYEK                    600318754.4  559494606 572010358.4
      YSGVCHTDLHAWHGDWPLPVK            104098406.5  118865215 115105378.7
      ANGTVVLVGLPAGAK                    4592718.9    4567092   5291864.9
      CSSDVFNHVVK                               NA    3767568   3911322.3
      DIPVPKPKPNELLINVK                         NA         NA          NA
      VVGLSSLPEIYEK                     36047939.9   33445820  32690717.2
      DIPVPEPKPNEILINVK                   769151.7   19057414  20505064.8
      EALDFFSR                           7209705.9         NA   7059163.8
      GVIFYENK                                  NA    5506509   6046983.4
      IQQGTDLAEVAPILCAGVTVYK             9431467.9    6314130   9385626.2
      IVGLSELPK                         18116936.8   17865395  19302307.0
      NMVSDIQEATK                        5769852.4    5392012   6477628.9
      VLGIDAGEEK                                NA         NA          NA
                                       50000amol_R1 50000amol_R2 50000amol_R3
      AAADALSDLEIK                               NA           NA           NA
      AAADALSDLEIKDSK                       4079331      2612596     10016845
      AEWALR                                4639152      4837731      4848298
      DEGLHTDFACLLFAHLK                          NA      4936813      2375920
      DIHDWNNR                              2208952      1448557           NA
      ELETLREENR                            4886849      4954543      5460366
      ESEFLFNAIHTIPEIGEK                   35045508     34834450     30415325
      GMMPGLTFSNELICR                       8409502      7324606      5882992
      IVTEAVEIEQR                          20994411     22811599     21164020
      LLVAFGNK                             14530874     15182336     13665955
      LLVAFGNKK                                  NA           NA           NA
      NKPDPAIVEK                           17949196     21424083     20858293
      TNFFEK                                5147855      5730161      5412857
      TVLFPIK                              14079305     14220523     13671086
      VENPFDFMENISLAGK                     11137935     12390987     12505961
      WIQDADALFGER                         12539068     12165644     15252099
      YFLDALPVALLGMNADLMNQYVEFVADR         22227514     23006993     22951679
      AANLGGVAVSGLEMAQNSQK                  9262947      8722671      8418522
      DAVWFGPPK                                  NA           NA           NA
      EIGYLFGAYR                                 NA           NA           NA
      FHPSVNLSILK                          12132653      5794715      7464436
      FLGFEQIFK                            37333151     39830328     38467650
      GANIASFVMVADAMLDQGDVF                41350390     40115097     43197863
      GCIISETGITSEQIHDIASAK                 7160191      7424782      9087030
      GGLCVDLK                             10795035     12816723     11135646
      ICYAFMR                               8574581      8566609      8504472
      NSWEGVLTGK                           13809137     14032908     13095732
      SLEEIVDEYSTFSESK                      5263726      6793873      5418863
      VLPIVSVPER                           42676778     41781054     37445935
      VTISGSGNVAQYAALK                      4807770      5937384           NA
      VTWENDNGEQEVAQGYR                          NA      2784590           NA
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK           NA           NA           NA
      AANLGGVAVSGLEMAQNSQR                 27841658     22170468     23410381
      ALVAQGVK                             14135946     14185583     13756692
      FIAEGSNMGSTPEAIAVFETAR               18507032     29152110     32916127
      GANIASFIK                            16688372     17875979     17635927
      GCIISETGITSEQVADISSAK                 3367072      4115433      4560523
      HIGQDTDVPAGDIGVGGR                   40588502     44567505     36003038
      IMINCFNECIDYAK                       13536503     10754948     10364695
      ITWTSER                              13939888     16325990           NA
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR        73349517      9756539      9492125
      SLEQIVNEYSTFSENK                     52865702     46265005     54447208
      STATGPSEAVWYGPPK                     39672610           NA     40401059
      VDIALPCATQNEVSGEEAK                  52868671     49810620     44882296
      VIELGGTVVSLSDSK                      12247965     12155890     12508867
      VQYIAGARPWTHVQK                       4032185      4292706           NA
      VTWENDKGEQEVAQGYR                    32832580     38567562     39569708
      AAGLTAAYAR                           63260717     62985078     66208504
      APEAEQVLSAAATFPIAQPATDVEAR           24253889     28456790     25120472
      AVQDNGESAFR                          12513677     12501178     14390444
      DGKAPEAEQVLSAAATFPIAQPATDVEAR        12757598     11345896     13845863
      GFTLAEVK                             29186051     28816679     32256844
      IAPRPLDLLRPVVR                       16813479     10327970     11056077
      IIVFPR                               55184418     55057804     50625874
      NQEIFDANVQR                          79784994     82495726     85988298
      TIGIAVDHR                                  NA     23402668     25824317
      VHFDQAGK                             17312631     21908353     21209596
      VHFDQAGKK                                  NA           NA           NA
      ANELLINVK                           220475704    211911441    211704216
      ANGTTVLVGMPAGAK                      18556953     17578887     18587105
      ATDGGAHGVINVSVSEAAIEASTR             97877743    117164178    150103877
      CCSDVFNQVVK                          48818841     51891192     50027633
      DIVGAVLK                             54771746     53777047     55469668
      EALDFFAR                            194442429    192641975    191456556
      EKDIVGAVLK                           39576094     43044430     39747974
      GVIFYESHGK                          109252083    116312387    133986343
      IGDYAGIK                            146129580    147510235    152473501
      LPLVGGHEGAGVVVGMGENVK               136107477    139805495    117906023
      SANLMAGHWVAISGAAGGLGSLAVQYAK         18029537     18218885     34675563
      SIGGEVFIDFTK                         54033367     58821698     45468016
      SIPETQK                              11754713      9457002      8723290
      SISIVGSYVGNR                        303039706    305804253    295256416
      VLGIDGGEGK                           19756621     23296503     25672498
      VLGIDGGEGKEELFR                     329061174    336926861    317481597
      VVGLSTLPEIYEK                       520312113    483464290    530396398
      YSGVCHTDLHAWHGDWPLPVK               102427714     89773629    105330304
      ANGTVVLVGLPAGAK                       5564601      6057040      3788527
      CSSDVFNHVVK                           3920848      3956395      3778704
      DIPVPKPKPNELLINVK                          NA           NA           NA
      VVGLSSLPEIYEK                        33256348     35422869     36327325
      DIPVPEPKPNEILINVK                    17797020     18331858     17240958
      EALDFFSR                              7001667      6931978           NA
      GVIFYENK                                   NA      6749750      4441632
      IQQGTDLAEVAPILCAGVTVYK                6661163      6580942      5993306
      IVGLSELPK                            18316412     17729357     18501551
      NMVSDIQEATK                           5167094      4095025      5952219
      VLGIDAGEEK                                 NA           NA           NA
                                       5000amol_R1 5000amol_R2 5000amol_R3 500amol_R1
      AAADALSDLEIK                              NA          NA          NA         NA
      AAADALSDLEIKDSK                           NA     2151018          NA   12292213
      AEWALR                               5389184     4814858     5324903    5283206
      DEGLHTDFACLLFAHLK                    4905130          NA          NA         NA
      DIHDWNNR                                  NA          NA          NA    1651494
      ELETLREENR                                NA     4388923          NA    3975954
      ESEFLFNAIHTIPEIGEK                  27620085    31077102    33854202   30725140
      GMMPGLTFSNELICR                     11530261     9728905     8894435   10021118
      IVTEAVEIEQR                         15187980    15704112    17652421   18159542
      LLVAFGNK                            10264512          NA    11037065    7157726
      LLVAFGNKK                                 NA     1471194          NA    1948976
      NKPDPAIVEK                          14621397    15485741    17175432   17885995
      TNFFEK                               6043369     7232939     6299788    6539659
      TVLFPIK                             11193640    11620405    10335729    9211933
      VENPFDFMENISLAGK                    17540354    15612370    17777556   16864253
      WIQDADALFGER                        17896679    19711235    19458892   16223464
      YFLDALPVALLGMNADLMNQYVEFVADR        21726912    21571079    23041208   26159856
      AANLGGVAVSGLEMAQNSQK                 8840594     5688446     9655673   10151091
      DAVWFGPPK                          170315067          NA          NA         NA
      EIGYLFGAYR                                NA          NA          NA         NA
      FHPSVNLSILK                          5163727     4674190     5118042    4519674
      FLGFEQIFK                           40224871    42117191    38411864   40934131
      GANIASFVMVADAMLDQGDVF               53891058    49400286    53615911   68059811
      GCIISETGITSEQIHDIASAK                7189590     7439801     6838449    6409398
      GGLCVDLK                            10737061          NA     8721433   10026126
      ICYAFMR                              9445433     9588446     8864904   10994691
      NSWEGVLTGK                          10048026    11635599    12147646   12897296
      SLEEIVDEYSTFSESK                     8685997     8018919     8776664    9603249
      VLPIVSVPER                          42250433    43435478    42450176   39023641
      VTISGSGNVAQYAALK                          NA     6885032     3292450         NA
      VTWENDNGEQEVAQGYR                    3289081     3378964     3837936    2802008
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK          NA          NA          NA         NA
      AANLGGVAVSGLEMAQNSQR                33928484    29018783     8754068   32505277
      ALVAQGVK                            11230058    12577107    11913989   12215333
      FIAEGSNMGSTPEAIAVFETAR              39045599    37380083    37558918   37712079
      GANIASFIK                           17675052    15443931    18320012   19640838
      GCIISETGITSEQVADISSAK                3954780     4168380     2857425    3629038
      HIGQDTDVPAGDIGVGGR                  31382936    29924490    27443385   33168669
      IMINCFNECIDYAK                      14110532    11150756    12917424   12478652
      ITWTSER                             16320367    19910662    17634413   16977371
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR       81401060    90920609    79238917   91432192
      SLEQIVNEYSTFSENK                    61810433    56717886    54699451   52707656
      STATGPSEAVWYGPPK                    43201405    46557318    49608204   41887488
      VDIALPCATQNEVSGEEAK                 51802294    47526864    54858568   55522308
      VIELGGTVVSLSDSK                     12766989    13577256    12713626   14120432
      VQYIAGARPWTHVQK                     12977317     4181010     9817608    5847700
      VTWENDKGEQEVAQGYR                   33277305    34213458    34900109   38003162
      AAGLTAAYAR                          56033262    58318432    56615755   62736795
      APEAEQVLSAAATFPIAQPATDVEAR          26243240    25674072    27229856   18513349
      AVQDNGESAFR                         13611358    13628645    11966821   12062708
      DGKAPEAEQVLSAAATFPIAQPATDVEAR       30024194    26252241    26770185   22882580
      GFTLAEVK                            26823470    28691011    28303193   29271013
      IAPRPLDLLRPVVR                      17715949    15757748     9308218   28623812
      IIVFPR                              52232962    51778690    54936664   46365093
      NQEIFDANVQR                         98320784    88211833    94988678   92932288
      TIGIAVDHR                           26190513    26312231    24548898   21780610
      VHFDQAGK                             8955267     8921302     9072787   10886933
      VHFDQAGKK                                 NA          NA          NA         NA
      ANELLINVK                          171758073   171725927   168410559  177339024
      ANGTTVLVGMPAGAK                     20578610    17057463    19884103         NA
      ATDGGAHGVINVSVSEAAIEASTR            85290609    89457390    98871172   81922308
      CCSDVFNQVVK                         48593705    48923833    49748445   46799733
      DIVGAVLK                                  NA    52726674          NA   56710896
      EALDFFAR                           207700885   207804703   201040422  210208493
      EKDIVGAVLK                          30155401    32558428    34653560   35142412
      GVIFYESHGK                          91251046    87568458    96484187   96437363
      IGDYAGIK                           143765894   148323874   135937546  126599559
      LPLVGGHEGAGVVVGMGENVK              130244431   130698486   133180797  142430928
      SANLMAGHWVAISGAAGGLGSLAVQYAK        23650998    13868940    22650372   37219371
      SIGGEVFIDFTK                        65797992    63540203    59340090   54177141
      SIPETQK                             32407296    31408177    17799399   23310612
      SISIVGSYVGNR                       288088494   287007312   285994834  284583764
      VLGIDGGEGK                          17847223    19439369    19969402    9702741
      VLGIDGGEGKEELFR                    319734613   321288225   319690018  340154662
      VVGLSTLPEIYEK                      586918292   575997334   585532550  560488576
      YSGVCHTDLHAWHGDWPLPVK              126864698   137441141   153218286  123342177
      ANGTVVLVGLPAGAK                      3821546     4780765     4413496    4543271
      CSSDVFNHVVK                               NA     3607608     1417890    3957113
      DIPVPKPKPNELLINVK                         NA          NA          NA         NA
      VVGLSSLPEIYEK                       38300000    39091900    38272265   33791673
      DIPVPEPKPNEILINVK                   18147523    19777858    20048843         NA
      EALDFFSR                             6407367     6504126     6727176    6340061
      GVIFYENK                                  NA          NA          NA    4880920
      IQQGTDLAEVAPILCAGVTVYK               7170633     8694551     8183542    8377114
      IVGLSELPK                           17028165    18767657    18354347   18373494
      NMVSDIQEATK                          5098747     5690121     7334190    5817290
      VLGIDAGEEK                                NA          NA     4619178         NA
                                       500amol_R2  500amol_R3 50amol_R1   50amol_R2
      AAADALSDLEIK                             NA          NA  17115971          NA
      AAADALSDLEIKDSK                     3690983   3227057.8   1636373          NA
      AEWALR                              5495591   5051812.5   3876689   4537838.0
      DEGLHTDFACLLFAHLK                        NA          NA   8667341   6082289.0
      DIHDWNNR                                 NA          NA        NA   1328446.1
      ELETLREENR                          3602883   6036040.8   5554150   6172061.8
      ESEFLFNAIHTIPEIGEK                 31897919  30811054.6  29899374  32126419.8
      GMMPGLTFSNELICR                     8086528   7588087.4   4898780          NA
      IVTEAVEIEQR                        16416107  16109719.9  26246064  22820544.0
      LLVAFGNK                            8461783   8604333.4   7065878   6061881.0
      LLVAFGNKK                                NA          NA   1044768   1664161.5
      NKPDPAIVEK                         19045469  18867421.8  19850979  16405210.7
      TNFFEK                              5506894   6112976.3   4401999   5601589.4
      TVLFPIK                             9287207  11415890.7  10083158   9580481.5
      VENPFDFMENISLAGK                   14194182  15709946.9  13255588  13889354.8
      WIQDADALFGER                       14500665  13826230.8  15354062  13350663.2
      YFLDALPVALLGMNADLMNQYVEFVADR       23711975  22037320.5  27007157  28163120.7
      AANLGGVAVSGLEMAQNSQK                9318267   9902346.0        NA   9792286.0
      DAVWFGPPK                                NA 170025230.1        NA          NA
      EIGYLFGAYR                               NA          NA  26549562  22735272.4
      FHPSVNLSILK                         5259949   4713024.1   9962350   6551595.5
      FLGFEQIFK                          44313121  41978593.5  51254289  55263513.9
      GANIASFVMVADAMLDQGDVF              68544250  68593193.1 108177056  93551788.6
      GCIISETGITSEQIHDIASAK               8897987   7992424.0   8523724   9227799.8
      GGLCVDLK                            9557098          NA  11662489   9969680.0
      ICYAFMR                            12277417  10837116.7   9808210   9353994.7
      NSWEGVLTGK                         11723159  13736622.9  12224789  13579369.5
      SLEEIVDEYSTFSESK                    8026155   8207733.6   9001774   6963744.5
      VLPIVSVPER                         40693598  38695941.9  28264686  30135305.0
      VTISGSGNVAQYAALK                         NA    131796.9   1906361   1948557.3
      VTWENDNGEQEVAQGYR                   2839483   2921613.7   1941576   2079875.5
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK         NA          NA   6139778   4457961.2
      AANLGGVAVSGLEMAQNSQR                8743736  34147046.6  26252876  26722784.3
      ALVAQGVK                           13613950  12981751.6        NA    675716.3
      FIAEGSNMGSTPEAIAVFETAR             38893108  39146423.4  34318717  38572422.9
      GANIASFIK                          20041844  20284595.7  15551729  17179239.3
      GCIISETGITSEQVADISSAK               3684328   3630532.8   4833739   3859823.6
      HIGQDTDVPAGDIGVGGR                 38255284  31550316.0  37466470  35010985.1
      IMINCFNECIDYAK                     12871923  12649747.8  13231756  13842815.7
      ITWTSER                            18388294  18250099.5        NA  18310895.6
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR      96818577  98114418.7  83341616          NA
      SLEQIVNEYSTFSENK                   48393955  48559277.2  42523736   9554806.2
      STATGPSEAVWYGPPK                   44969165  42091348.7  36189475  37469954.0
      VDIALPCATQNEVSGEEAK                51672579  51176822.3   7099913  44241661.1
      VIELGGTVVSLSDSK                    14610908  12391960.6  19458506  14962765.1
      VQYIAGARPWTHVQK                     9095806          NA        NA          NA
      VTWENDKGEQEVAQGYR                  37092438  36823624.9  50429801  48489701.1
      AAGLTAAYAR                         65268886  64198038.9  45087204  58283660.7
      APEAEQVLSAAATFPIAQPATDVEAR         16998084  12929383.6   4555803   2690555.8
      AVQDNGESAFR                        13582235  13883834.5  17192983  16439196.0
      DGKAPEAEQVLSAAATFPIAQPATDVEAR      26500427  24508826.6  27289193  22241342.8
      GFTLAEVK                           31580998  32359325.2  27287393  33072821.0
      IAPRPLDLLRPVVR                     22765877  25874254.6  18867693   9867749.7
      IIVFPR                             46262994  46150456.3   1739506   1942643.1
      NQEIFDANVQR                        90494706  86754854.0  68225808  82192945.6
      TIGIAVDHR                          21016977  22591181.4  29146164  28911252.0
      VHFDQAGK                           10690211   9881791.7  13339836  14228700.4
      VHFDQAGKK                                NA          NA   1273995   1253359.1
      ANELLINVK                         166607876 149129934.1 162766276 174134329.2
      ANGTTVLVGMPAGAK                          NA          NA  19323839  18408395.4
      ATDGGAHGVINVSVSEAAIEASTR           94835045  97315513.3 104913796  96564392.4
      CCSDVFNQVVK                        52721555  49518761.7  57121323  61367647.9
      DIVGAVLK                           46201062  47480177.2  47797478  56351914.0
      EALDFFAR                          214059869 206520679.3 202335466 206202866.2
      EKDIVGAVLK                         33660848  31767521.6  48436642  31242679.5
      GVIFYESHGK                        100219648  95558566.7 194495272 217144614.9
      IGDYAGIK                          125000149 132680132.9        NA 125070387.1
      LPLVGGHEGAGVVVGMGENVK             120715728 147068476.8 109544427 105604089.7
      SANLMAGHWVAISGAAGGLGSLAVQYAK       27546541  34045891.5  48905979  40595723.4
      SIGGEVFIDFTK                       67877632  75527570.3  73335090  41837687.2
      SIPETQK                            19088656  18048143.6  27575573  28468039.9
      SISIVGSYVGNR                      286546486 278116805.3 274724203 264159058.4
      VLGIDGGEGK                          8126008   8854655.4  10315412   8624609.9
      VLGIDGGEGKEELFR                   364237972 366524890.0 423773507 462057731.5
      VVGLSTLPEIYEK                     599301040 567932840.6 545624528 535957682.6
      YSGVCHTDLHAWHGDWPLPVK              90394614 115450972.0  72221250 118848979.3
      ANGTVVLVGLPAGAK                     5849371   4972980.4   4151931   4312172.9
      CSSDVFNHVVK                              NA   1104099.2   3884026   3066825.8
      DIPVPKPKPNELLINVK                        NA          NA   5194680   4603664.9
      VVGLSSLPEIYEK                      30247315  32676223.7  26549686  28177753.5
      DIPVPEPKPNEILINVK                  18649233  19881822.6   3779971  25262383.3
      EALDFFSR                            6609484   9770902.6   5771613   5433546.1
      GVIFYENK                            5093531   6019715.9   5020517   5708752.7
      IQQGTDLAEVAPILCAGVTVYK              8231946   7581328.5   5979160   6845802.9
      IVGLSELPK                          20807835  19619218.1  14475083  15733894.2
      NMVSDIQEATK                         6717867   4740682.3   5325652   5886637.1
      VLGIDAGEEK                               NA          NA        NA          NA
                                       50amol_R3
      AAADALSDLEIK                            NA
      AAADALSDLEIKDSK                         NA
      AEWALR                             4622274
      DEGLHTDFACLLFAHLK                  4960864
      DIHDWNNR                           1424338
      ELETLREENR                         6299854
      ESEFLFNAIHTIPEIGEK                26142175
      GMMPGLTFSNELICR                         NA
      IVTEAVEIEQR                       20824333
      LLVAFGNK                           7202535
      LLVAFGNKK                          2304657
      NKPDPAIVEK                        20294259
      TNFFEK                             5111491
      TVLFPIK                           11001068
      VENPFDFMENISLAGK                  14664253
      WIQDADALFGER                      14662454
      YFLDALPVALLGMNADLMNQYVEFVADR      26459835
      AANLGGVAVSGLEMAQNSQK               8054903
      DAVWFGPPK                               NA
      EIGYLFGAYR                        25626429
      FHPSVNLSILK                        7147875
      FLGFEQIFK                         51642133
      GANIASFVMVADAMLDQGDVF             85540881
      GCIISETGITSEQIHDIASAK              7369392
      GGLCVDLK                          11255623
      ICYAFMR                            9768861
      NSWEGVLTGK                        14686249
      SLEEIVDEYSTFSESK                   5699260
      VLPIVSVPER                        32028573
      VTISGSGNVAQYAALK                        NA
      VTWENDNGEQEVAQGYR                  2505580
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK   4111844
      AANLGGVAVSGLEMAQNSQR              27956956
      ALVAQGVK                                NA
      FIAEGSNMGSTPEAIAVFETAR            37556128
      GANIASFIK                         17208058
      GCIISETGITSEQVADISSAK                   NA
      HIGQDTDVPAGDIGVGGR                31981282
      IMINCFNECIDYAK                    13378997
      ITWTSER                           17225686
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     98073519
      SLEQIVNEYSTFSENK                  43899292
      STATGPSEAVWYGPPK                  36602916
      VDIALPCATQNEVSGEEAK               44489200
      VIELGGTVVSLSDSK                   13800126
      VQYIAGARPWTHVQK                         NA
      VTWENDKGEQEVAQGYR                 46591640
      AAGLTAAYAR                        54084543
      APEAEQVLSAAATFPIAQPATDVEAR         2835897
      AVQDNGESAFR                       14564614
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     28538039
      GFTLAEVK                          32227580
      IAPRPLDLLRPVVR                     9717348
      IIVFPR                            40616550
      NQEIFDANVQR                       79538102
      TIGIAVDHR                         24478797
      VHFDQAGK                          14976680
      VHFDQAGKK                          1176479
      ANELLINVK                        172164178
      ANGTTVLVGMPAGAK                   18750410
      ATDGGAHGVINVSVSEAAIEASTR          84758926
      CCSDVFNQVVK                       62377182
      DIVGAVLK                          53032619
      EALDFFAR                         200765180
      EKDIVGAVLK                        29761029
      GVIFYESHGK                       159904734
      IGDYAGIK                         120474497
      LPLVGGHEGAGVVVGMGENVK            127556185
      SANLMAGHWVAISGAAGGLGSLAVQYAK      39166263
      SIGGEVFIDFTK                      59986676
      SIPETQK                           52989140
      SISIVGSYVGNR                     259374198
      VLGIDGGEGK                         9185334
      VLGIDGGEGKEELFR                  468313620
      VVGLSTLPEIYEK                    538799355
      YSGVCHTDLHAWHGDWPLPVK             99179130
      ANGTVVLVGLPAGAK                    4504047
      CSSDVFNHVVK                        4152785
      DIPVPKPKPNELLINVK                  3427055
      VVGLSSLPEIYEK                     31438774
      DIPVPEPKPNEILINVK                 19354857
      EALDFFSR                           5224860
      GVIFYENK                           5557418
      IQQGTDLAEVAPILCAGVTVYK             7373382
      IVGLSELPK                         15482663
      NMVSDIQEATK                        6232538
      VLGIDAGEEK                              NA

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
                                       125amol_R1  125amol_R2  125amol_R3
      AAADALSDLEIK                             NA          NA          NA
      AAADALSDLEIKDSK                    13218450  14668357.4  13788629.4
      AEWALR                              6207258   6424893.7   5476021.2
      DEGLHTDFACLLFAHLK                   3077795   5325002.9   4292601.1
      DIHDWNNR                            2126158          NA   2625754.2
      ELETLREENR                          5318855   5325422.5   6854984.7
      ESEFLFNAIHTIPEIGEK                 34560611  34920614.3  41209976.4
      GMMPGLTFSNELICR                     6759773   8298117.3   8229321.2
      IVTEAVEIEQR                        18004130  20813479.2  19350150.5
      LLVAFGNK                            6662884   8153982.8   7955419.5
      LLVAFGNKK                           1706084   2228629.4   1877478.6
      NKPDPAIVEK                         19822461  20341423.1  17749792.4
      TNFFEK                              5933181   5660477.5   5130802.5
      TVLFPIK                            11876033  11427954.3  10808201.2
      VENPFDFMENISLAGK                   19185382  18986097.5  15865775.5
      WIQDADALFGER                       17080744  18002122.7  17703521.8
      YFLDALPVALLGMNADLMNQYVEFVADR       25939182  24538526.5  21947337.6
      AANLGGVAVSGLEMAQNSQK               10161986  10379989.7  10110119.2
      DAVWFGPPK                                NA          NA          NA
      EIGYLFGAYR                         25282194  28075800.5  24609907.1
      FHPSVNLSILK                         5020603   5295315.8   5266194.2
      FLGFEQIFK                          50978679  49844931.0  48901955.5
      GANIASFVMVADAMLDQGDVF              82452488  91526438.1  82439064.6
      GCIISETGITSEQIHDIASAK               8556773   8338714.1   8986951.4
      GGLCVDLK                           11219045          NA          NA
      ICYAFMR                            12005914  14651573.2  12307971.2
      NSWEGVLTGK                         13496224  19407800.9          NA
      SLEEIVDEYSTFSESK                    9566620  10805889.3   9368784.2
      VLPIVSVPER                         34406081  42320356.4  33621603.2
      VTISGSGNVAQYAALK                    6961989    125220.8          NA
      VTWENDNGEQEVAQGYR                        NA   3270090.1   2436849.5
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK         NA          NA   1911678.6
      AANLGGVAVSGLEMAQNSQR                9922133  38491456.7  35513667.1
      ALVAQGVK                           14857602  15426794.3  14331805.6
      FIAEGSNMGSTPEAIAVFETAR             42101280  46106246.5  44989074.8
      GANIASFIK                          19170213  20602627.4  20223256.0
      GCIISETGITSEQVADISSAK               3007167   3399643.3   2827735.3
      HIGQDTDVPAGDIGVGGR                 35657487  38815601.9  34521867.5
      IMINCFNECIDYAK                     12909391  14141752.6  12862212.1
      ITWTSER                            17884678  18682932.5    550911.5
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR       9749591 122619201.7 102992249.9
      SLEQIVNEYSTFSENK                   47646336  56047748.6  46676944.7
      STATGPSEAVWYGPPK                   39445835  47229740.1  48589126.2
      VDIALPCATQNEVSGEEAK                48477763  57223692.8  51690260.2
      VIELGGTVVSLSDSK                    13512341  14419741.2  13641770.6
      VQYIAGARPWTHVQK                     3799819          NA   3900206.6
      VTWENDKGEQEVAQGYR                  44446718  45035203.6  46200156.6
      AAGLTAAYAR                         60414470  65885398.3  63966043.1
      APEAEQVLSAAATFPIAQPATDVEAR         20326816  26831670.3  21973490.6
      AVQDNGESAFR                        14902160  14849836.8  15155622.9
      DGKAPEAEQVLSAAATFPIAQPATDVEAR      41710690  47044064.7  42789211.1
      GFTLAEVK                           32933783  36923181.3  33192091.6
      IAPRPLDLLRPVVR                     23824116  30963734.9  28780293.5
      IIVFPR                             42209356  48123499.7   3359746.2
      NQEIFDANVQR                        81425766  86907631.2  82135288.3
      TIGIAVDHR                          20481345  21278192.2  18073686.4
      VHFDQAGK                                 NA          NA  14186958.6
      VHFDQAGKK                                NA          NA          NA
      ANELLINVK                         164758979 188958822.9 171593397.0
      ANGTTVLVGMPAGAK                          NA          NA          NA
      ATDGGAHGVINVSVSEAAIEASTR           82960635  68544647.8  83175370.2
      CCSDVFNQVVK                        57960007  64652807.3  56210693.5
      DIVGAVLK                                 NA  62858994.0  52919447.9
      EALDFFAR                          198110845 218173852.1 203892265.5
      EKDIVGAVLK                         32617139  36758486.2  35871761.1
      GVIFYESHGK                        117783851 126143887.5 111190187.6
      IGDYAGIK                          135550025 144323293.5 141507468.7
      LPLVGGHEGAGVVVGMGENVK             113185882  90758560.1 152089346.5
      SANLMAGHWVAISGAAGGLGSLAVQYAK       23233490  31786161.6  33891179.6
      SIGGEVFIDFTK                       81313898  67515565.4  78098684.1
      SIPETQK                            27784058  51081718.1  49200702.4
      SISIVGSYVGNR                      252148833 297888396.9 279041703.8
      VLGIDGGEGK                         10133545  13010915.9  11214577.5
      VLGIDGGEGKEELFR                   403929214 412503624.8 331719740.0
      VVGLSTLPEIYEK                     528036458 579758355.9 545308707.2
      YSGVCHTDLHAWHGDWPLPVK              91661699 112831904.7 119126595.9
      ANGTVVLVGLPAGAK                     4787860   6578993.8   5169931.3
      CSSDVFNHVVK                         4185669   5506482.3   4351847.5
      DIPVPKPKPNELLINVK                        NA   1378403.9          NA
      VVGLSSLPEIYEK                      36980944  39750273.0  36808237.1
      DIPVPEPKPNEILINVK                  18493316  23719246.9  21442384.9
      EALDFFSR                            6922266   8048766.8   6661553.6
      GVIFYENK                                 NA   6644976.8   5866705.7
      IQQGTDLAEVAPILCAGVTVYK              8018858   8335672.0   8397806.4
      IVGLSELPK                          19105747  21319103.7  21847755.3
      NMVSDIQEATK                         3736205   4790215.8          NA
      VLGIDAGEEK                               NA          NA          NA
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
                                       2500amol_R1 2500amol_R2 2500amol_R3
      AAADALSDLEIK                              NA          NA          NA
      AAADALSDLEIKDSK                    6252233.7   8474606.1     8946931
      AEWALR                             6070218.4   5284894.3     5016018
      DEGLHTDFACLLFAHLK                  5367006.8   5702891.6     3942675
      DIHDWNNR                            683491.6          NA          NA
      ELETLREENR                                NA          NA     4779228
      ESEFLFNAIHTIPEIGEK                37374713.7  29614298.8    37005204
      GMMPGLTFSNELICR                    9086000.5  10881900.9     9664109
      IVTEAVEIEQR                       19100979.5  20655069.9    19176343
      LLVAFGNK                          11694000.5  11933706.8     9894192
      LLVAFGNKK                           828201.1          NA          NA
      NKPDPAIVEK                        15649699.0  16584792.5    18762606
      TNFFEK                             7227937.2   6496855.5     6357103
      TVLFPIK                           11057902.9  13110334.5    12727826
      VENPFDFMENISLAGK                  18140992.6  17923076.6    16833210
      WIQDADALFGER                      20834053.3  24613342.4    19891355
      YFLDALPVALLGMNADLMNQYVEFVADR      26872201.4  26712566.3    25659941
      AANLGGVAVSGLEMAQNSQK              11685505.0  11851723.9     8703641
      DAVWFGPPK                                 NA          NA          NA
      EIGYLFGAYR                                NA          NA          NA
      FHPSVNLSILK                        5905300.3   5376923.0     8064877
      FLGFEQIFK                         43136759.3  47888400.8    40930030
      GANIASFVMVADAMLDQGDVF             54503684.6  63443211.1    59678192
      GCIISETGITSEQIHDIASAK              6609892.5   9126774.3     8887192
      GGLCVDLK                          10586086.3  10719667.2    11707412
      ICYAFMR                            9182317.8  10009302.8    10120768
      NSWEGVLTGK                        13793441.0  15213022.5    12261469
      SLEEIVDEYSTFSESK                   9104584.3  12397891.6     8619346
      VLPIVSVPER                        44483289.7  44767277.5    44648503
      VTISGSGNVAQYAALK                   7442341.4    197867.4          NA
      VTWENDNGEQEVAQGYR                  3658888.4   3252757.8     3307109
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK          NA          NA          NA
      AANLGGVAVSGLEMAQNSQR              34578645.4  38515793.3    34187460
      ALVAQGVK                                  NA  14162255.8    13942211
      FIAEGSNMGSTPEAIAVFETAR            39169379.5  42332039.2    37671723
      GANIASFIK                         19435488.3  19993433.3    19152613
      GCIISETGITSEQVADISSAK              4039060.2   3606900.8     3725180
      HIGQDTDVPAGDIGVGGR                36930825.9  32213502.9    33445623
      IMINCFNECIDYAK                    14073791.2  15574440.0    14164040
      ITWTSER                                   NA  19247504.5    17386235
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     95615338.6 107817894.4    87880396
      SLEQIVNEYSTFSENK                  56757105.3  59326746.2    61366155
      STATGPSEAVWYGPPK                  52933086.4  56182529.1    46753299
      VDIALPCATQNEVSGEEAK               57721340.0  63880838.0    56617983
      VIELGGTVVSLSDSK                   14962628.7  16125226.5    14205310
      VQYIAGARPWTHVQK                    8894534.1          NA     6569956
      VTWENDKGEQEVAQGYR                 36795960.4  36379619.4    35243985
      AAGLTAAYAR                        59032826.6  62290831.8    58583491
      APEAEQVLSAAATFPIAQPATDVEAR        28631823.5  31747008.7    30101679
      AVQDNGESAFR                       14357327.1  14791561.1    14136182
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     30450913.9  25616766.8    21977579
      GFTLAEVK                          32414427.0  35864628.3    31290270
      IAPRPLDLLRPVVR                     9479977.4  22366318.4    10578664
      IIVFPR                            63949574.0  63787308.3    53628966
      NQEIFDANVQR                       92314852.5  97582734.5    90831237
      TIGIAVDHR                         24075071.3  26064786.0    24250553
      VHFDQAGK                           9906661.9          NA    11276136
      VHFDQAGKK                                 NA          NA          NA
      ANELLINVK                        165523415.7 163746392.0   160263714
      ANGTTVLVGMPAGAK                           NA          NA          NA
      ATDGGAHGVINVSVSEAAIEASTR          76760727.4 103595198.1    93652077
      CCSDVFNQVVK                       56174104.4  58016174.6    55658444
      DIVGAVLK                          51081076.2          NA    51319881
      EALDFFAR                         208722829.0 221388447.4   205124760
      EKDIVGAVLK                        31887708.5  34634884.9    37846091
      GVIFYESHGK                        91446191.7  96966130.8    95441153
      IGDYAGIK                         146971455.3 144024310.5   148584159
      LPLVGGHEGAGVVVGMGENVK            108625073.4 129151921.2   159799421
      SANLMAGHWVAISGAAGGLGSLAVQYAK      40473432.6  24060246.6    29850961
      SIGGEVFIDFTK                      56072158.8  79989898.9    76279255
      SIPETQK                           43250386.1  24239223.3    24836594
      SISIVGSYVGNR                     280955476.4 293764329.6   266896399
      VLGIDGGEGK                        20904140.9  22034922.8    22030199
      VLGIDGGEGKEELFR                  297999492.7 342850142.5   345011159
      VVGLSTLPEIYEK                    565075511.5 579053244.8   523846604
      YSGVCHTDLHAWHGDWPLPVK            157527064.2 169566022.6   148480983
      ANGTVVLVGLPAGAK                    5357555.6   4455366.0     5318221
      CSSDVFNHVVK                               NA          NA          NA
      DIPVPKPKPNELLINVK                         NA          NA          NA
      VVGLSSLPEIYEK                     39342474.5  43207293.1    36082808
      DIPVPEPKPNEILINVK                 23606759.0  26726422.6    23038231
      EALDFFSR                           7725665.0   7233546.8     9221999
      GVIFYENK                           5358829.9   5848266.9     5918191
      IQQGTDLAEVAPILCAGVTVYK             7199265.1   9527797.6     7275166
      IVGLSELPK                         20685382.8  22046469.7    20501127
      NMVSDIQEATK                        6761536.5   7681681.5     6387540
      VLGIDAGEEK                                NA          NA          NA
                                        250amol_R1 250amol_R2  250amol_R3
      AAADALSDLEIK                              NA         NA          NA
      AAADALSDLEIKDSK                   11025078.7   14410611  13471469.8
      AEWALR                             5324150.3    6410555   5476389.6
      DEGLHTDFACLLFAHLK                         NA         NA          NA
      DIHDWNNR                                  NA    1913884   2340844.1
      ELETLREENR                         4275483.1    4197643   4417215.1
      ESEFLFNAIHTIPEIGEK                28762650.8   25587626  39910244.1
      GMMPGLTFSNELICR                    7870974.3         NA   8248658.3
      IVTEAVEIEQR                       18815451.9   17958059  20834207.0
      LLVAFGNK                           7233202.2    6734732   7848744.1
      LLVAFGNKK                                 NA    1824942          NA
      NKPDPAIVEK                        17118137.8   20709119  20293724.7
      TNFFEK                             6934786.2    5984043   5880576.4
      TVLFPIK                            9706461.8   10308044  12522955.3
      VENPFDFMENISLAGK                  18249325.2   18813740  19610110.6
      WIQDADALFGER                      15067927.4   15650684  18498169.6
      YFLDALPVALLGMNADLMNQYVEFVADR      24846674.2   24307683  26203781.9
      AANLGGVAVSGLEMAQNSQK              11148752.3   10898945   6426506.4
      DAVWFGPPK                                 NA         NA          NA
      EIGYLFGAYR                                NA   21794822  23459712.3
      FHPSVNLSILK                        5010381.7    6876732   5806233.4
      FLGFEQIFK                         48504588.9   47222967  47792100.7
      GANIASFVMVADAMLDQGDVF             77710759.0   87501320  83419075.4
      GCIISETGITSEQIHDIASAK              9002481.6    9345224   7555858.2
      GGLCVDLK                          11874801.9   12075634  11516760.8
      ICYAFMR                           10753209.9   11375345  13125390.2
      NSWEGVLTGK                        12737318.9   14667213  15331113.7
      SLEEIVDEYSTFSESK                   9057068.5    8684704   8836246.0
      VLPIVSVPER                        34545459.7   34806901  42749378.8
      VTISGSGNVAQYAALK                   3441432.3         NA    172954.4
      VTWENDNGEQEVAQGYR                  3046956.0    2816489   3471370.9
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK          NA         NA          NA
      AANLGGVAVSGLEMAQNSQR              38011520.1   36195210  42635087.1
      ALVAQGVK                          14383457.7   13675564  13102959.1
      FIAEGSNMGSTPEAIAVFETAR            44085388.5   28660718  46149290.7
      GANIASFIK                         18353808.1   19058074  20685734.6
      GCIISETGITSEQVADISSAK              4310239.6    3834515   3801214.7
      HIGQDTDVPAGDIGVGGR                34413256.8   38250072  33986300.9
      IMINCFNECIDYAK                    15028479.7   14008840  14206355.0
      ITWTSER                           18113923.8   17841514  17935256.1
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR    108054513.3  102487497 105462424.8
      SLEQIVNEYSTFSENK                  55879162.8   55793069  64365469.4
      STATGPSEAVWYGPPK                  48125038.8   47051217  41072388.0
      VDIALPCATQNEVSGEEAK               50379950.4   44890035  54465883.1
      VIELGGTVVSLSDSK                   14201145.7   13891273  14223445.3
      VQYIAGARPWTHVQK                   11026144.8   10574869   4935693.6
      VTWENDKGEQEVAQGYR                 37723659.1   43867716  42238804.5
      AAGLTAAYAR                        64746349.7   59801582  63270618.8
      APEAEQVLSAAATFPIAQPATDVEAR        20840074.7   19263560  16115133.6
      AVQDNGESAFR                         286453.7   14537379  12610543.4
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     28312734.6   28258946  31391771.9
      GFTLAEVK                          31734016.7   35858867  33736354.5
      IAPRPLDLLRPVVR                    13439913.1   25521176  26020060.7
      IIVFPR                            44904193.3   50141688  47092464.5
      NQEIFDANVQR                       93233933.3   85388186  83934990.4
      TIGIAVDHR                         21808140.8   24116509  22360589.1
      VHFDQAGK                          11227647.6   11446908  11164272.3
      VHFDQAGKK                                 NA         NA          NA
      ANELLINVK                        168953120.2  175695772 170006274.8
      ANGTTVLVGMPAGAK                           NA         NA          NA
      ATDGGAHGVINVSVSEAAIEASTR          79985927.5   79958649  70722012.2
      CCSDVFNQVVK                       53333191.7   59148320  59895274.3
      DIVGAVLK                          55180833.2   57496252          NA
      EALDFFAR                         196470505.4  190222928 218350607.6
      EKDIVGAVLK                        35114784.9   38251094  40084352.1
      GVIFYESHGK                       102289829.0   94612572 107487631.5
      IGDYAGIK                         138215955.7  134608765          NA
      LPLVGGHEGAGVVVGMGENVK            122276982.1  126440435  91813642.6
      SANLMAGHWVAISGAAGGLGSLAVQYAK      27732747.9   31829907  31803008.5
      SIGGEVFIDFTK                      65806104.9   63070959  57491943.6
      SIPETQK                           28859670.7   24754436  27119183.9
      SISIVGSYVGNR                     252368858.9  261304779 277664809.9
      VLGIDGGEGK                        10425048.3   11393747  11463353.5
      VLGIDGGEGKEELFR                  325549537.7  338949920 370241109.5
      VVGLSTLPEIYEK                    548545906.8  509329633 565284749.7
      YSGVCHTDLHAWHGDWPLPVK            101948660.4  116810188 117528213.2
      ANGTVVLVGLPAGAK                    4902807.0    4872782   5797901.9
      CSSDVFNHVVK                               NA    4033151   4433557.7
      DIPVPKPKPNELLINVK                         NA         NA          NA
      VVGLSSLPEIYEK                     38145855.3   35375310  36621419.5
      DIPVPEPKPNEILINVK                   786553.7   20506700  22737645.0
      EALDFFSR                           7646229.5         NA   7494760.2
      GVIFYENK                                  NA    5865250   6537059.6
      IQQGTDLAEVAPILCAGVTVYK             9972573.3    6734323   9931203.4
      IVGLSELPK                         18471084.8   19161328  21528502.6
      NMVSDIQEATK                        6141995.5    5737971   6963570.7
      VLGIDAGEEK                                NA         NA          NA
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
                                       5000amol_R1 5000amol_R2 5000amol_R3 500amol_R1
      AAADALSDLEIK                              NA          NA          NA         NA
      AAADALSDLEIKDSK                           NA     2289399          NA   12894341
      AEWALR                               5848626     5048690     5464988    5679358
      DEGLHTDFACLLFAHLK                    5425231          NA          NA         NA
      DIHDWNNR                                  NA          NA          NA    1961258
      ELETLREENR                                NA     4564929          NA    4314729
      ESEFLFNAIHTIPEIGEK                  29330027    33466564    36643775   32026500
      GMMPGLTFSNELICR                     12035649    10330899     9409919   10443837
      IVTEAVEIEQR                         16531288    17190632    19350869   18909041
      LLVAFGNK                            10605354          NA    11917068    7453387
      LLVAFGNKK                                 NA     1432927          NA    2311791
      NKPDPAIVEK                          15989934    16959153    18852490   18670732
      TNFFEK                               6480772     7664505     6519640    6932239
      TVLFPIK                             11768967    12452960    11462726    9635864
      VENPFDFMENISLAGK                    19181129    17095797    19595763   17678567
      WIQDADALFGER                        19614612    21591792    21423870   16966488
      YFLDALPVALLGMNADLMNQYVEFVADR        23447042    23486453    25254087   27221386
      AANLGGVAVSGLEMAQNSQK                 9160777     6026408    10319354   10571062
      DAVWFGPPK                          163784706          NA          NA         NA
      EIGYLFGAYR                                NA          NA          NA         NA
      FHPSVNLSILK                          5541093     4950898     5294100    4794576
      FLGFEQIFK                           41582418    44610188    41301045   42422865
      GANIASFVMVADAMLDQGDVF               54327221    51644296    57489783   69597729
      GCIISETGITSEQIHDIASAK                7514540     7904958     7223494    6679783
      GGLCVDLK                            11174675          NA     9590044   10469472
      ICYAFMR                              9761061    10162031     9620441   11454041
      NSWEGVLTGK                          10805116    12618156    13498134   13548505
      SLEEIVDEYSTFSESK                     9045714     8510678     9307343   10004246
      VLPIVSVPER                          44117395    46242778    45688717   40553418
      VTISGSGNVAQYAALK                          NA     6648853     3188876         NA
      VTWENDNGEQEVAQGYR                    3792678     3381531     3802506    3170749
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK          NA          NA          NA         NA
      AANLGGVAVSGLEMAQNSQR                36434556    31527035     9562333   33857020
      ALVAQGVK                            11861856    13512460    13228536   12793700
      FIAEGSNMGSTPEAIAVFETAR              40970148    39914326    40433177   39223253
      GANIASFIK                           19303982    16911225    20097364   20490808
      GCIISETGITSEQVADISSAK                4431716     4241369     2870837    4026289
      HIGQDTDVPAGDIGVGGR                  33182434    32107262    29604166   34549162
      IMINCFNECIDYAK                      15083015    12071580    14334852   13110813
      ITWTSER                             17838729    21789620    19494798   17827629
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR       81868381    94918636    84837279   93360306
      SLEQIVNEYSTFSENK                    63129725    59669581    58912957   54290862
      STATGPSEAVWYGPPK                    44609808    49293815    53337337   43400788
      VDIALPCATQNEVSGEEAK                 53132645    50111639    59057960   57339512
      VIELGGTVVSLSDSK                     13851484    14800381    14139367   14854934
      VQYIAGARPWTHVQK                     13949367     4420026    10147714    6215507
      VTWENDKGEQEVAQGYR                   34821479    36467633    37547088   39500490
      AAGLTAAYAR                          56606503    61041120    60781880   64237192
      APEAEQVLSAAATFPIAQPATDVEAR          28562942    28103797    29849060   19275525
      AVQDNGESAFR                         14342898    14622948    13279019   12625649
      DGKAPEAEQVLSAAATFPIAQPATDVEAR       32357417    28567774    29317383   23816695
      GFTLAEVK                            28590909    30994735    30715853   30520689
      IAPRPLDLLRPVVR                      19422840    17262013    10246423   29924441
      IIVFPR                              54981441    55402006    59185777   48257172
      NQEIFDANVQR                         97198102    90684711   100089407   93888194
      TIGIAVDHR                           28342205    28688103    26944711   22655530
      VHFDQAGK                             9468810     9590370    10076179   11405619
      VHFDQAGKK                                 NA          NA          NA         NA
      ANELLINVK                          164923348   170733624   173852781  173823422
      ANGTTVLVGMPAGAK                     22404285    18677400    21804099         NA
      ATDGGAHGVINVSVSEAAIEASTR            83873990    91492338   103718769   82479726
      CCSDVFNQVVK                         49478003    51403639    53575786   48121402
      DIVGAVLK                                  NA    55391808          NA   58306042
      EALDFFAR                           198553248   205506464   206891899  204832123
      EKDIVGAVLK                          31750142    34842182    37336566   36578114
      GVIFYESHGK                          89052815    88828800   100630750   96472569
      IGDYAGIK                           138994265   148727821   140931814  125240592
      LPLVGGHEGAGVVVGMGENVK              126479188   131790083   138450658  141703877
      SANLMAGHWVAISGAAGGLGSLAVQYAK        25289845    15024722    24650447   38790309
      SIGGEVFIDFTK                        66306937    66413065    63614264   55389363
      SIPETQK                             35019244    34220146    19522724   24252488
      SISIVGSYVGNR                       272435152   280283366   291756030  273533580
      VLGIDGGEGK                          19490760    21272616    22087660   10191286
      VLGIDGGEGKEELFR                    300182067   311160832   323989577  324091245
      VVGLSTLPEIYEK                      541903842   546932923   583544283  523378640
      YSGVCHTDLHAWHGDWPLPVK              123652559   139213727   159664006  123227780
      ANGTVVLVGLPAGAK                      4169628     4985161     4512374    4914300
      CSSDVFNHVVK                               NA     3670115     1422422    4394861
      DIPVPKPKPNELLINVK                         NA          NA          NA         NA
      VVGLSSLPEIYEK                       40292955    41813065    41225859   35167247
      DIPVPEPKPNEILINVK                   19898274    21666233    22068326         NA
      EALDFFSR                             6800391     6923262     7012112    6652724
      GVIFYENK                                  NA          NA          NA    5217361
      IQQGTDLAEVAPILCAGVTVYK               7534716     9245904     8595541    8740922
      IVGLSELPK                           18578847    20549627    20127439   19158744
      NMVSDIQEATK                          5482663     6017639     7575689    6185695
      VLGIDAGEEK                                NA          NA     4698300         NA
                                       500amol_R2  500amol_R3 50amol_R1   50amol_R2
      AAADALSDLEIK                             NA          NA  16458044          NA
      AAADALSDLEIKDSK                     3729110   3625605.2   1745328          NA
      AEWALR                              5376759   5481807.8   4528760   5588844.9
      DEGLHTDFACLLFAHLK                        NA          NA  10922386   7824804.0
      DIHDWNNR                                 NA          NA        NA   2057646.2
      ELETLREENR                          3482088   6463946.5   6697556   7744526.6
      ESEFLFNAIHTIPEIGEK                 34877920  36028414.5  30216421  33442814.7
      GMMPGLTFSNELICR                     8333672   8701601.4   4847661          NA
      IVTEAVEIEQR                        17532258  19286424.1  24127591  23349098.6
      LLVAFGNK                            8867548  10041364.8   6404343   6552082.2
      LLVAFGNKK                                NA          NA   1861951   2565306.9
      NKPDPAIVEK                         20216142  22674392.8  17866190  16701727.0
      TNFFEK                              5483882   6762109.0   4906051   6720051.4
      TVLFPIK                             9836395  13510389.5   8375589   9882267.1
      VENPFDFMENISLAGK                   15007899  18915183.5  11595333  14064564.0
      WIQDADALFGER                       15350152  16660083.5  13592110  13549214.2
      YFLDALPVALLGMNADLMNQYVEFVADR       25581009  26007047.7  26047144  29115890.2
      AANLGGVAVSGLEMAQNSQK                9698141  11470162.2        NA  10810445.9
      DAVWFGPPK                                NA 184662909.9        NA          NA
      EIGYLFGAYR                               NA          NA  25474273  23486403.1
      FHPSVNLSILK                         5233187   5208323.8  11130416   7870297.5
      FLGFEQIFK                          47894177  48749519.0  51796427  56498996.8
      GANIASFVMVADAMLDQGDVF              69278233  77506299.6 104251616  90894705.1
      GCIISETGITSEQIHDIASAK               9156046   9151344.9   8493788  10467184.5
      GGLCVDLK                           10089701          NA   9923997  10393040.0
      ICYAFMR                            12897338  12683874.1   8727822  10010418.0
      NSWEGVLTGK                         12440033  16383514.6  10244711  13756086.5
      SLEEIVDEYSTFSESK                    8296159   9440678.5   8798742   7825261.6
      VLPIVSVPER                         44389465  45073280.4  28712928  31126944.4
      VTISGSGNVAQYAALK                         NA    122200.9   3522570   3055850.3
      VTWENDNGEQEVAQGYR                   2568543   2929715.4   2827717   2893007.1
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK         NA          NA   7534011   5658321.0
      AANLGGVAVSGLEMAQNSQR                9476256  40110996.6  25750858  27722706.6
      ALVAQGVK                           14422975  15395312.8        NA    692362.8
      FIAEGSNMGSTPEAIAVFETAR             42583451  45643437.8  34872678  39989496.2
      GANIASFIK                          21303811  24358271.7  14052298  17509993.2
      GCIISETGITSEQVADISSAK               3443028   3755805.3   6411041   5107546.7
      HIGQDTDVPAGDIGVGGR                 41913447  36864503.5  38005742  36458938.2
      IMINCFNECIDYAK                     13658223  15075131.1  11061861  14055410.4
      ITWTSER                            19443655  21965951.0        NA  18510485.8
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR      97281119 110414804.4  79863200          NA
      SLEQIVNEYSTFSENK                   50973489  56036045.1  42403537   9553651.5
      STATGPSEAVWYGPPK                   48537110  48859082.8  36527465  38242066.8
      VDIALPCATQNEVSGEEAK                55048622  59263391.1   7098566  44605675.6
      VIELGGTVVSLSDSK                    15471971  14852812.3  16465136  15127301.4
      VQYIAGARPWTHVQK                     9024172          NA        NA          NA
      VTWENDKGEQEVAQGYR                  40538634  42916043.2  51260591  50184817.7
      AAGLTAAYAR                         66351679  72809947.3  43710031  56878872.8
      APEAEQVLSAAATFPIAQPATDVEAR         18156822  15476156.0   4192286   2753505.7
      AVQDNGESAFR                        14387203  16447338.2  14286333  16903107.0
      DGKAPEAEQVLSAAATFPIAQPATDVEAR      28624549  28894001.3  26448233  23005837.2
      GFTLAEVK                           34461222  37877171.1  27382794  34436899.7
      IAPRPLDLLRPVVR                     24103713  31181013.1  16731478  10017740.9
      IIVFPR                             50683444  53842641.2   1766762   2017919.4
      NQEIFDANVQR                        89319381  95945997.4  64031611  78713044.3
      TIGIAVDHR                          22627383  26700597.4  27867411  29863742.4
      VHFDQAGK                           11326260  11724392.8  11088653  14557945.1
      VHFDQAGKK                                NA          NA   2235618   1919792.4
      ANELLINVK                         164387881 161750242.4 145723745 163740267.8
      ANGTTVLVGMPAGAK                          NA          NA  17870918  18842762.8
      ATDGGAHGVINVSVSEAAIEASTR           93512466 107329997.8  98144662  92488353.4
      CCSDVFNQVVK                        54906500  56882771.4  56587431  60969633.4
      DIVGAVLK                           48053702  54513852.6  47312740  55936963.5
      EALDFFAR                          210708080 222914512.7 178937695 192601686.9
      EKDIVGAVLK                         36877298  37065547.8  49192106  32457883.3
      GVIFYESHGK                         98944044 104995118.2 180600363 207741803.3
      IGDYAGIK                          123597875 144858266.6        NA 118750139.4
      LPLVGGHEGAGVVVGMGENVK             119334213 161112006.5 100886882 100735780.8
      SANLMAGHWVAISGAAGGLGSLAVQYAK       29971325  39908888.6  48631842  42175384.8
      SIGGEVFIDFTK                       68537652  85288522.2  70603481  40619266.5
      SIPETQK                            20573864  21314950.4  26484481  29417960.2
      SISIVGSYVGNR                      280376908 296513745.1 236035746 243038229.6
      VLGIDGGEGK                          8591425  10657475.2   8885421   8714354.4
      VLGIDGGEGKEELFR                   355382564 387653826.4 357702031 421579929.3
      VVGLSTLPEIYEK                     584421493 590634105.6 443995263 482350928.7
      YSGVCHTDLHAWHGDWPLPVK              89276841 126764279.8  66942461 113642409.7
      ANGTVVLVGLPAGAK                     5676279   5348735.7   4956719   5377670.5
      CSSDVFNHVVK                              NA   1140953.0   5165616   4065954.2
      DIPVPKPKPNELLINVK                        NA          NA   6754685   6036640.9
      VVGLSSLPEIYEK                      33132816  38115445.8  26969885  29260517.7
      DIPVPEPKPNEILINVK                  19746269  23960438.1   3354018  25649407.7
      EALDFFSR                            6674131  10971274.5   6168418   6381000.7
      GVIFYENK                            5019810   6584041.2   5752751   6954384.6
      IQQGTDLAEVAPILCAGVTVYK              8416335   8624694.0   6102147   7862700.0
      IVGLSELPK                          22145908  23547712.2  13120101  16046326.6
      NMVSDIQEATK                         6663882   5222045.8   5998369   7102490.4
      VLGIDAGEEK                               NA          NA        NA          NA
                                       50amol_R3
      AAADALSDLEIK                            NA
      AAADALSDLEIKDSK                         NA
      AEWALR                             5232325
      DEGLHTDFACLLFAHLK                  5850027
      DIHDWNNR                           2023635
      ELETLREENR                         7260462
      ESEFLFNAIHTIPEIGEK                27362483
      GMMPGLTFSNELICR                         NA
      IVTEAVEIEQR                       21236910
      LLVAFGNK                           7262483
      LLVAFGNKK                          3262991
      NKPDPAIVEK                        20444913
      TNFFEK                             5637327
      TVLFPIK                           10769137
      VENPFDFMENISLAGK                  14572927
      WIQDADALFGER                      14637242
      YFLDALPVALLGMNADLMNQYVEFVADR      27366158
      AANLGGVAVSGLEMAQNSQK               8253766
      DAVWFGPPK                               NA
      EIGYLFGAYR                        26445530
      FHPSVNLSILK                        7894335
      FLGFEQIFK                         53646306
      GANIASFVMVADAMLDQGDVF             85397848
      GCIISETGITSEQIHDIASAK              7709291
      GGLCVDLK                          11064068
      ICYAFMR                            9775006
      NSWEGVLTGK                        14259620
      SLEEIVDEYSTFSESK                   5921693
      VLPIVSVPER                        33444872
      VTISGSGNVAQYAALK                        NA
      VTWENDNGEQEVAQGYR                  3181770
      YVAGARPWTHVSNVDIALPCATQNEVSGDEAK   4798197
      AANLGGVAVSGLEMAQNSQR              29069594
      ALVAQGVK                                NA
      FIAEGSNMGSTPEAIAVFETAR            39298483
      GANIASFIK                         17376155
      GCIISETGITSEQVADISSAK                   NA
      HIGQDTDVPAGDIGVGGR                33477948
      IMINCFNECIDYAK                    12993527
      ITWTSER                           16984751
      SEPEFQQAYEEVVSSLEDSTLFEQHPEYR     97391732
      SLEQIVNEYSTFSENK                  44893913
      STATGPSEAVWYGPPK                  37987369
      VDIALPCATQNEVSGEEAK               45836592
      VIELGGTVVSLSDSK                   13455678
      VQYIAGARPWTHVQK                         NA
      VTWENDKGEQEVAQGYR                 48708729
      AAGLTAAYAR                        54219631
      APEAEQVLSAAATFPIAQPATDVEAR         2893729
      AVQDNGESAFR                       14230218
      DGKAPEAEQVLSAAATFPIAQPATDVEAR     29570335
      GFTLAEVK                          33692945
      IAPRPLDLLRPVVR                     9707016
      IIVFPR                            42535379
      NQEIFDANVQR                       78044772
      TIGIAVDHR                         25228134
      VHFDQAGK                          14589465
      VHFDQAGKK                          1662734
      ANELLINVK                        167973174
      ANGTTVLVGMPAGAK                   19157689
      ATDGGAHGVINVSVSEAAIEASTR          83228585
      CCSDVFNQVVK                       63463576
      DIVGAVLK                          53940319
      EALDFFAR                         195086664
      EKDIVGAVLK                        31166273
      GVIFYESHGK                       157195768
      IGDYAGIK                         118156384
      LPLVGGHEGAGVVVGMGENVK            125341321
      SANLMAGHWVAISGAAGGLGSLAVQYAK      40876964
      SIGGEVFIDFTK                      59844463
      SIPETQK                           54712102
      SISIVGSYVGNR                     249626034
      VLGIDGGEGK                         9044017
      VLGIDGGEGKEELFR                  448332731
      VVGLSTLPEIYEK                    511931369
      YSGVCHTDLHAWHGDWPLPVK             97501987
      ANGTVVLVGLPAGAK                    5159005
      CSSDVFNHVVK                        5037909
      DIPVPKPKPNELLINVK                  4124999
      VVGLSSLPEIYEK                     32925756
      DIPVPEPKPNEILINVK                 19338690
      EALDFFSR                           5635306
      GVIFYENK                           6225169
      IQQGTDLAEVAPILCAGVTVYK             7791614
      IVGLSELPK                         15663531
      NMVSDIQEATK                        6914262
      VLGIDAGEEK                              NA

