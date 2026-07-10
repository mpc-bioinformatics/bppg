# generate graphs from fasta

    Code
      igraph::as_edgelist(res[[1]])
    Output
           [,1]                  
      [1,] "sp|P09938|RIR2_YEAST"
           [,2]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
      [1,] "AAADALSDLEIK;ELETLR;DAENHK;SHQVHR;EEPLLNEDK;TVLFPIK;YHEIWQAYK;AEASFWTAEEIDLSK;DIHDWNNR;MNENER;VLAFFAASDGIVNENLVENFSTEVQIPEAK;SFYGFQIMIENIHSETYSLLIDTYIK;ESEFLFNAIHTIPEIGEK;AEWALR;WIQDADALFGER;LVAFASIEGVFFSGSFASIFWLK;GMMPGLTFSNELICR;DEGLHTDFACLLFAHLK;NKPDPAIVEK;IVTEAVEIEQR;YFLDALPVALLGMNADLMNQYVEFVADR;LLVAFGNK;VENPFDFMENISLAGK;TNFFEK;VSDYQK;AGVMSK;QEAGAFTFNEDF;MPKETPSK;ETPSKAAADALSDLEIK;AAADALSDLEIKDSK;DSKSNLNK;SNLNKELETLR;ELETLREENR;EENRVK;VKSDMLK;SDMLKEK;LSKDAENHK;DAENHKAYLK;AYLKSHQVHR;SHQVHRHK;LKEMEK;EMEKEEPLLNEDK;EEPLLNEDKER;ERTVLFPIK;TVLFPIKYHEIWQAYK;YHEIWQAYKR;RAEASFWTAEEIDLSK;AEASFWTAEEIDLSKDIHDWNNR;DIHDWNNRMNENER;MNENERFFISR;FFISRVLAFFAASDGIVNENLVENFSTEVQIPEAK;SFYGFQIMIENIHSETYSLLIDTYIKDPK;DPKESEFLFNAIHTIPEIGEK;ESEFLFNAIHTIPEIGEKAEWALR;AEWALRWIQDADALFGER;WIQDADALFGERLVAFASIEGVFFSGSFASIFWLK;LVAFASIEGVFFSGSFASIFWLKK;RGMMPGLTFSNELICR;GMMPGLTFSNELICRDEGLHTDFACLLFAHLK;DEGLHTDFACLLFAHLKNKPDPAIVEK;NKPDPAIVEKIVTEAVEIEQR;IVTEAVEIEQRYFLDALPVALLGMNADLMNQYVEFVADR;YFLDALPVALLGMNADLMNQYVEFVADRLLVAFGNK;LLVAFGNKK;YYKVENPFDFMENISLAGK;VENPFDFMENISLAGKTNFFEK;TNFFEKR;RVSDYQK;VSDYQKAGVMSK;AGVMSKSTK;STKQEAGAFTFNEDF;MPKETPSKAAADALSDLEIK;ETPSKAAADALSDLEIKDSK;AAADALSDLEIKDSKSNLNK;DSKSNLNKELETLR;SNLNKELETLREENR;ELETLREENRVK;EENRVKSDMLK;VKSDMLKEK;SDMLKEKLSK;EKLSKDAENHK;LSKDAENHKAYLK;DAENHKAYLKSHQVHR;AYLKSHQVHRHK;SHQVHRHKLK;HKLKEMEK;LKEMEKEEPLLNEDK;EMEKEEPLLNEDKER;EEPLLNEDKERTVLFPIK;ERTVLFPIKYHEIWQAYK;TVLFPIKYHEIWQAYKR;YHEIWQAYKRAEASFWTAEEIDLSK;RAEASFWTAEEIDLSKDIHDWNNR;AEASFWTAEEIDLSKDIHDWNNRMNENER;DIHDWNNRMNENERFFISR;MNENERFFISRVLAFFAASDGIVNENLVENFSTEVQIPEAK;SFYGFQIMIENIHSETYSLLIDTYIKDPKESEFLFNAIHTIPEIGEK;DPKESEFLFNAIHTIPEIGEKAEWALR;ESEFLFNAIHTIPEIGEKAEWALRWIQDADALFGER;AEWALRWIQDADALFGERLVAFASIEGVFFSGSFASIFWLK;WIQDADALFGERLVAFASIEGVFFSGSFASIFWLKK;LVAFASIEGVFFSGSFASIFWLKKR;KRGMMPGLTFSNELICR;RGMMPGLTFSNELICRDEGLHTDFACLLFAHLK;GMMPGLTFSNELICRDEGLHTDFACLLFAHLKNKPDPAIVEK;DEGLHTDFACLLFAHLKNKPDPAIVEKIVTEAVEIEQR;NKPDPAIVEKIVTEAVEIEQRYFLDALPVALLGMNADLMNQYVEFVADR;IVTEAVEIEQRYFLDALPVALLGMNADLMNQYVEFVADRLLVAFGNK;YFLDALPVALLGMNADLMNQYVEFVADRLLVAFGNKK;LLVAFGNKKYYK;KYYKVENPFDFMENISLAGK;YYKVENPFDFMENISLAGKTNFFEK;VENPFDFMENISLAGKTNFFEKR;TNFFEKRVSDYQK;RVSDYQKAGVMSK;VSDYQKAGVMSKSTK;AGVMSKSTKQEAGAFTFNEDF;PKETPSK;PKETPSKAAADALSDLEIK"

---

    Code
      igraph::as_edgelist(res[[2]])
    Output
           [,1]                   
      [1,] "sp|P40212|RL13B_YEAST"
      [2,] "sp|P40212|RL13B_YEAST"
      [3,] "sp|Q12690|RL13A_YEAST"
      [4,] "sp|Q12690|RL13A_YEAST"
           [,2]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
      [1,] "NARAAR;AARAAK;IIVFPRDGK;DGKAPEAEQVLSAAATFPIAQPATDVEAR;RNARAAR;NARAARAAK;AARAAKIAPRPLDLLRPVVR;EYQSKIIVFPRDGK;IIVFPRDGKAPEAEQVLSAAATFPIAQPATDVEAR;DGKAPEAEQVLSAAATFPIAQPATDVEARAVQDNGESAFR"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
      [2,] "NLPILK;VHFDQAGK;IAPRPLDLLRPVVR;GFTLAEVK;AAGLTAAYAR;TIGIAVDHR;NQEIFDANVQR;IIVFPR;APEAEQVLSAAATFPIAQPATDVEAR;AVQDNGESAFR;AEAEAEK;MAISKNLPILK;NLPILKNHFR;KHWQER;HWQERVK;VKVHFDQAGK;VHFDQAGKK;AAKIAPRPLDLLRPVVR;IAPRPLDLLRPVVRAPTVK;APTVKYNR;AGRGFTLAEVK;GFTLAEVKAAGLTAAYAR;AAGLTAAYARTIGIAVDHR;TIGIAVDHRR;QNRNQEIFDANVQR;NQEIFDANVQRLK;LKEYQSK;EYQSKIIVFPR;APEAEQVLSAAATFPIAQPATDVEARAVQDNGESAFR;AVQDNGESAFRTLR;TLRLAR;LARSEK;EKAEAEAEK;AEAEAEKK;MAISKNLPILKNHFR;NLPILKNHFRK;NHFRKHWQER;KHWQERVK;HWQERVKVHFDQAGK;VKVHFDQAGKK;VHFDQAGKKVSR;VSRRNAR;AAKIAPRPLDLLRPVVRAPTVK;IAPRPLDLLRPVVRAPTVKYNR;APTVKYNRK;YNRKVR;KVRAGR;VRAGRGFTLAEVK;AGRGFTLAEVKAAGLTAAYAR;GFTLAEVKAAGLTAAYARTIGIAVDHR;AAGLTAAYARTIGIAVDHRR;TIGIAVDHRRQNR;RQNRNQEIFDANVQR;QNRNQEIFDANVQRLK;NQEIFDANVQRLKEYQSK;LKEYQSKIIVFPR;APEAEQVLSAAATFPIAQPATDVEARAVQDNGESAFRTLR;AVQDNGESAFRTLRLAR;TLRLARSEK;LARSEKK;SEKKFR;KFRGIR;FRGIREK;GIREKR;AREKAEAEAEK;EKAEAEAEKK;AEAEAEKKK;AISKNLPILK;AISKNLPILKNHFR"
      [3,] "NLPILK;VHFDQAGK;IAPRPLDLLRPVVR;GFTLAEVK;AAGLTAAYAR;TIGIAVDHR;NQEIFDANVQR;IIVFPR;APEAEQVLSAAATFPIAQPATDVEAR;AVQDNGESAFR;AEAEAEK;MAISKNLPILK;NLPILKNHFR;KHWQER;HWQERVK;VKVHFDQAGK;VHFDQAGKK;AAKIAPRPLDLLRPVVR;IAPRPLDLLRPVVRAPTVK;APTVKYNR;AGRGFTLAEVK;GFTLAEVKAAGLTAAYAR;AAGLTAAYARTIGIAVDHR;TIGIAVDHRR;QNRNQEIFDANVQR;NQEIFDANVQRLK;LKEYQSK;EYQSKIIVFPR;APEAEQVLSAAATFPIAQPATDVEARAVQDNGESAFR;AVQDNGESAFRTLR;TLRLAR;LARSEK;EKAEAEAEK;AEAEAEKK;MAISKNLPILKNHFR;NLPILKNHFRK;NHFRKHWQER;KHWQERVK;HWQERVKVHFDQAGK;VKVHFDQAGKK;VHFDQAGKKVSR;VSRRNAR;AAKIAPRPLDLLRPVVRAPTVK;IAPRPLDLLRPVVRAPTVKYNR;APTVKYNRK;YNRKVR;KVRAGR;VRAGRGFTLAEVK;AGRGFTLAEVKAAGLTAAYAR;GFTLAEVKAAGLTAAYARTIGIAVDHR;AAGLTAAYARTIGIAVDHRR;TIGIAVDHRRQNR;RQNRNQEIFDANVQR;QNRNQEIFDANVQRLK;NQEIFDANVQRLKEYQSK;LKEYQSKIIVFPR;APEAEQVLSAAATFPIAQPATDVEARAVQDNGESAFRTLR;AVQDNGESAFRTLRLAR;TLRLARSEK;LARSEKK;SEKKFR;KFRGIR;FRGIREK;GIREKR;AREKAEAEAEK;EKAEAEAEKK;AEAEAEKKK;AISKNLPILK;AISKNLPILKNHFR"
      [4,] "NARATR;ATRAAK;IIVFPRNGK;NGKAPEAEQVLSAAATFPIAQPATDVEAR;RNARATR;NARATRAAK;ATRAAKIAPRPLDLLRPVVR;EYQSKIIVFPRNGK;IIVFPRNGKAPEAEQVLSAAATFPIAQPATDVEAR;NGKAPEAEQVLSAAATFPIAQPATDVEARAVQDNGESAFR"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

---

    Code
      igraph::as_edgelist(res[[3]])
    Output
           [,1]                  
      [1,] "sp|P07262|DHE4_YEAST"
      [2,] "sp|P39708|DHE5_YEAST"
      [3,] "sp|P07262|DHE4_YEAST"
      [4,] "sp|P39708|DHE5_YEAST"
           [,2]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
      [1,] "MSEPEFQQAYEEVVSSLEDSTLFEQHPEYR;VTWENDK;GEQEVAQGYR;VQYNSAK;NSLTGLDMGGGK;SNNEIR;HIGQDTDVPAGDIGVGGR;GLNWGGSLIRPEATGYGLVYYTQAMIDYATNGK;VIELGGTVVSLSDSK;GCIISETGITSEQVADISSAK;SLEQIVNEYSTFSENK;VQYIAGARPWTHVQK;VDIALPCATQNEVSGEEAK;ALVAQGVK;FIAEGSNMGSTPEAIAVFETAR;STATGPSEAVWYGPPK;AANLGGVAVSGLEMAQNSQR;ITWTSER;IMINCFNECIDYAK;VLPSLVK;GANIASFIK;VSDAMFDQGDVF;MSEPEFQQAYEEVVSSLEDSTLFEQHPEYRK;IIQFRVTWENDK;VTWENDKGEQEVAQGYR;GEQEVAQGYRVQYNSAK;VQYNSAKGPYK;FLGFEQIFKNSLTGLDMGGGK;NSLTGLDMGGGKGGLCVDLK;GGLCVDLKGR;GRSNNEIR;SNNEIRR;ELSRHIGQDTDVPAGDIGVGGR;HIGQDTDVPAGDIGVGGREIGYLFGAYR;NSWEGVLTGKGLNWGGSLIRPEATGYGLVYYTQAMIDYATNGK;GLNWGGSLIRPEATGYGLVYYTQAMIDYATNGKESFEGK;VTISGSGNVAQYAALKVIELGGTVVSLSDSK;VIELGGTVVSLSDSKGCIISETGITSEQVADISSAK;GCIISETGITSEQVADISSAKVNFK;VNFKSLEQIVNEYSTFSENK;SLEQIVNEYSTFSENKVQYIAGARPWTHVQK;VQYIAGARPWTHVQKVDIALPCATQNEVSGEEAK;VDIALPCATQNEVSGEEAKALVAQGVK;ALVAQGVKFIAEGSNMGSTPEAIAVFETAR;FIAEGSNMGSTPEAIAVFETARSTATGPSEAVWYGPPK;STATGPSEAVWYGPPKAANLGGVAVSGLEMAQNSQR;AANLGGVAVSGLEMAQNSQRITWTSER;ITWTSERVDQELK;VDQELKR;RIMINCFNECIDYAK;IMINCFNECIDYAKK;YTKDGK;DGKVLPSLVK;VLPSLVKGANIASFIK;GANIASFIKVSDAMFDQGDVF;MSEPEFQQAYEEVVSSLEDSTLFEQHPEYRKVLPIVSVPER;VLPIVSVPERIIQFRVTWENDK;IIQFRVTWENDKGEQEVAQGYR;VTWENDKGEQEVAQGYRVQYNSAK;GEQEVAQGYRVQYNSAKGPYK;VQYNSAKGPYKGGLR;FHPSVNLSILKFLGFEQIFKNSLTGLDMGGGK;FLGFEQIFKNSLTGLDMGGGKGGLCVDLK;NSLTGLDMGGGKGGLCVDLKGR;GGLCVDLKGRSNNEIR;GRSNNEIRR;SNNEIRRICYAFMR;ICYAFMRELSRHIGQDTDVPAGDIGVGGR;ELSRHIGQDTDVPAGDIGVGGREIGYLFGAYR;HIGQDTDVPAGDIGVGGREIGYLFGAYRSYK;SYKNSWEGVLTGKGLNWGGSLIRPEATGYGLVYYTQAMIDYATNGK;NSWEGVLTGKGLNWGGSLIRPEATGYGLVYYTQAMIDYATNGKESFEGK;GLNWGGSLIRPEATGYGLVYYTQAMIDYATNGKESFEGKR;RVTISGSGNVAQYAALKVIELGGTVVSLSDSK;VIELGGTVVSLSDSKGCIISETGITSEQVADISSAKVNFK;GCIISETGITSEQVADISSAKVNFKSLEQIVNEYSTFSENK;VNFKSLEQIVNEYSTFSENKVQYIAGARPWTHVQK;SLEQIVNEYSTFSENKVQYIAGARPWTHVQKVDIALPCATQNEVSGEEAK;VQYIAGARPWTHVQKVDIALPCATQNEVSGEEAKALVAQGVK;VDIALPCATQNEVSGEEAKALVAQGVKFIAEGSNMGSTPEAIAVFETAR;ALVAQGVKFIAEGSNMGSTPEAIAVFETARSTATGPSEAVWYGPPK;STATGPSEAVWYGPPKAANLGGVAVSGLEMAQNSQRITWTSER;AANLGGVAVSGLEMAQNSQRITWTSERVDQELK;ITWTSERVDQELKR;VDQELKRIMINCFNECIDYAK;RIMINCFNECIDYAKK;IMINCFNECIDYAKKYTK;KYTKDGK;YTKDGKVLPSLVK;DGKVLPSLVKGANIASFIK;VLPSLVKGANIASFIKVSDAMFDQGDVF;SEPEFQQAYEEVVSSLEDSTLFEQHPEYR;SEPEFQQAYEEVVSSLEDSTLFEQHPEYRK;SEPEFQQAYEEVVSSLEDSTLFEQHPEYRKVLPIVSVPER"
      [2,] "VLPIVSVPER;FHPSVNLSILK;FLGFEQIFK;GGLCVDLK;ICYAFMR;EIGYLFGAYR;NSWEGVLTGK;ESFEGK;VTISGSGNVAQYAALK;VDQELK;KVLPIVSVPER;VLPIVSVPERIIQFR;GPYKGGLR;GGLRFHPSVNLSILK;FHPSVNLSILKFLGFEQIFK;RICYAFMR;ICYAFMRELSR;EIGYLFGAYRSYK;SYKNSWEGVLTGK;ESFEGKR;RVTISGSGNVAQYAALK;KVLPIVSVPERIIQFR;GPYKGGLRFHPSVNLSILK;GGLRFHPSVNLSILKFLGFEQIFK;RICYAFMRELSR;EIGYLFGAYRSYKNSWEGVLTGK;ESFEGKRVTISGSGNVAQYAALK"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
      [3,] "VLPIVSVPER;FHPSVNLSILK;FLGFEQIFK;GGLCVDLK;ICYAFMR;EIGYLFGAYR;NSWEGVLTGK;ESFEGK;VTISGSGNVAQYAALK;VDQELK;KVLPIVSVPER;VLPIVSVPERIIQFR;GPYKGGLR;GGLRFHPSVNLSILK;FHPSVNLSILKFLGFEQIFK;RICYAFMR;ICYAFMRELSR;EIGYLFGAYRSYK;SYKNSWEGVLTGK;ESFEGKR;RVTISGSGNVAQYAALK;KVLPIVSVPERIIQFR;GPYKGGLRFHPSVNLSILK;GGLRFHPSVNLSILKFLGFEQIFK;RICYAFMRELSR;EIGYLFGAYRSYKNSWEGVLTGK;ESFEGKRVTISGSGNVAQYAALK"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
      [4,] "MTSEPEFQQAYDEIVSSVEDSK;VTWENDNGEQEVAQGYR;VQFNSAK;NALTGLDMGGGK;SDNEIR;DTDVPAGDIGVGGR;GLNWGGSLIRPEATGFGLVYYTQAMIDYATNGK;VIELGGIVVSLSDSK;GCIISETGITSEQIHDIASAK;SLEEIVDEYSTFSESK;YVAGARPWTHVSNVDIALPCATQNEVSGDEAK;ALVASGVK;FVAEGANMGSTPEAISVFETAR;STATNAK;DAVWFGPPK;AANLGGVAVSGLEMAQNSQK;VTWTAER;IMINCFNDCIQAAQEYSTEK;NTNTLPSLVK;GANIASFVMVADAMLDQGDVF;MTSEPEFQQAYDEIVSSVEDSKIFEK;IFEKFPQYK;FPQYKK;IIQFRVTWENDNGEQEVAQGYR;VTWENDNGEQEVAQGYRVQFNSAK;VQFNSAKGPYK;FLGFEQIFKNALTGLDMGGGK;NALTGLDMGGGKGGLCVDLK;GGLCVDLKGK;GKSDNEIR;SDNEIRR;ELSRHIGK;HIGKDTDVPAGDIGVGGR;DTDVPAGDIGVGGREIGYLFGAYR;NSWEGVLTGKGLNWGGSLIRPEATGFGLVYYTQAMIDYATNGK;GLNWGGSLIRPEATGFGLVYYTQAMIDYATNGKESFEGK;VTISGSGNVAQYAALKVIELGGIVVSLSDSK;VIELGGIVVSLSDSKGCIISETGITSEQIHDIASAK;GCIISETGITSEQIHDIASAKIR;FKSLEEIVDEYSTFSESK;SLEEIVDEYSTFSESKMK;MKYVAGARPWTHVSNVDIALPCATQNEVSGDEAK;YVAGARPWTHVSNVDIALPCATQNEVSGDEAKALVASGVK;ALVASGVKFVAEGANMGSTPEAISVFETAR;FVAEGANMGSTPEAISVFETARSTATNAK;STATNAKDAVWFGPPK;DAVWFGPPKAANLGGVAVSGLEMAQNSQK;AANLGGVAVSGLEMAQNSQKVTWTAER;VTWTAERVDQELK;VDQELKK;KIMINCFNDCIQAAQEYSTEK;IMINCFNDCIQAAQEYSTEKNTNTLPSLVK;NTNTLPSLVKGANIASFVMVADAMLDQGDVF;MTSEPEFQQAYDEIVSSVEDSKIFEKFPQYK;IFEKFPQYKK;FPQYKKVLPIVSVPER;VLPIVSVPERIIQFRVTWENDNGEQEVAQGYR;IIQFRVTWENDNGEQEVAQGYRVQFNSAK;VTWENDNGEQEVAQGYRVQFNSAKGPYK;VQFNSAKGPYKGGLR;FHPSVNLSILKFLGFEQIFKNALTGLDMGGGK;FLGFEQIFKNALTGLDMGGGKGGLCVDLK;NALTGLDMGGGKGGLCVDLKGK;GGLCVDLKGKSDNEIR;GKSDNEIRR;SDNEIRRICYAFMR;ICYAFMRELSRHIGK;ELSRHIGKDTDVPAGDIGVGGR;HIGKDTDVPAGDIGVGGREIGYLFGAYR;DTDVPAGDIGVGGREIGYLFGAYRSYK;SYKNSWEGVLTGKGLNWGGSLIRPEATGFGLVYYTQAMIDYATNGK;NSWEGVLTGKGLNWGGSLIRPEATGFGLVYYTQAMIDYATNGKESFEGK;GLNWGGSLIRPEATGFGLVYYTQAMIDYATNGKESFEGKR;RVTISGSGNVAQYAALKVIELGGIVVSLSDSK;VIELGGIVVSLSDSKGCIISETGITSEQIHDIASAKIR;GCIISETGITSEQIHDIASAKIRFK;IRFKSLEEIVDEYSTFSESK;FKSLEEIVDEYSTFSESKMK;SLEEIVDEYSTFSESKMKYVAGARPWTHVSNVDIALPCATQNEVSGDEAK;MKYVAGARPWTHVSNVDIALPCATQNEVSGDEAKALVASGVK;ALVASGVKFVAEGANMGSTPEAISVFETARSTATNAK;FVAEGANMGSTPEAISVFETARSTATNAKDAVWFGPPK;STATNAKDAVWFGPPKAANLGGVAVSGLEMAQNSQK;DAVWFGPPKAANLGGVAVSGLEMAQNSQKVTWTAER;AANLGGVAVSGLEMAQNSQKVTWTAERVDQELK;VTWTAERVDQELKK;VDQELKKIMINCFNDCIQAAQEYSTEK;KIMINCFNDCIQAAQEYSTEKNTNTLPSLVK;TSEPEFQQAYDEIVSSVEDSK;TSEPEFQQAYDEIVSSVEDSKIFEK;TSEPEFQQAYDEIVSSVEDSKIFEKFPQYK"                                                                              

# test generateGraphsFromQuantData

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                          
       [1,] "sp|P09938|RIR2_YEAST" "AAADALSDLEIKDSK"             
       [2,] "sp|P09938|RIR2_YEAST" "AEWALR"                      
       [3,] "sp|P09938|RIR2_YEAST" "DEGLHTDFACLLFAHLK"           
       [4,] "sp|P09938|RIR2_YEAST" "DIHDWNNR"                    
       [5,] "sp|P09938|RIR2_YEAST" "ESEFLFNAIHTIPEIGEK"          
       [6,] "sp|P09938|RIR2_YEAST" "GMMPGLTFSNELICR"             
       [7,] "sp|P09938|RIR2_YEAST" "IVTEAVEIEQR"                 
       [8,] "sp|P09938|RIR2_YEAST" "LLVAFGNK"                    
       [9,] "sp|P09938|RIR2_YEAST" "NKPDPAIVEK"                  
      [10,] "sp|P09938|RIR2_YEAST" "TNFFEK"                      
      [11,] "sp|P09938|RIR2_YEAST" "TVLFPIK"                     
      [12,] "sp|P09938|RIR2_YEAST" "VENPFDFMENISLAGK"            
      [13,] "sp|P09938|RIR2_YEAST" "WIQDADALFGER"                
      [14,] "sp|P09938|RIR2_YEAST" "YFLDALPVALLGMNADLMNQYVEFVADR"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA  0.942555857  0.191337293  0.473638399  0.281388730
       [6]  0.004830101 -0.143762647 -0.037945653 -0.742235799  0.143241202
      [11] -0.272304407 -0.294216237 -0.070510314 -0.120199738 -0.304583227

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                    [,2]                           
       [1,] "sp|P40212|RL13B_YEAST" "AAGLTAAYAR"                   
       [2,] "sp|Q12690|RL13A_YEAST" "AAGLTAAYAR"                   
       [3,] "sp|P40212|RL13B_YEAST" "APEAEQVLSAAATFPIAQPATDVEAR"   
       [4,] "sp|Q12690|RL13A_YEAST" "APEAEQVLSAAATFPIAQPATDVEAR"   
       [5,] "sp|P40212|RL13B_YEAST" "AVQDNGESAFR"                  
       [6,] "sp|Q12690|RL13A_YEAST" "AVQDNGESAFR"                  
       [7,] "sp|P40212|RL13B_YEAST" "DGKAPEAEQVLSAAATFPIAQPATDVEAR"
       [8,] "sp|P40212|RL13B_YEAST" "GFTLAEVK"                     
       [9,] "sp|Q12690|RL13A_YEAST" "GFTLAEVK"                     
      [10,] "sp|P40212|RL13B_YEAST" "IAPRPLDLLRPVVR"               
      [11,] "sp|Q12690|RL13A_YEAST" "IAPRPLDLLRPVVR"               
      [12,] "sp|P40212|RL13B_YEAST" "IIVFPR"                       
      [13,] "sp|Q12690|RL13A_YEAST" "IIVFPR"                       
      [14,] "sp|P40212|RL13B_YEAST" "NQEIFDANVQR"                  
      [15,] "sp|Q12690|RL13A_YEAST" "NQEIFDANVQR"                  
      [16,] "sp|P40212|RL13B_YEAST" "TIGIAVDHR"                    
      [17,] "sp|Q12690|RL13A_YEAST" "TIGIAVDHR"                    
      [18,] "sp|P40212|RL13B_YEAST" "VHFDQAGK"                     
      [19,] "sp|Q12690|RL13A_YEAST" "VHFDQAGK"                     

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA  0.01301067 -0.39541482  0.25095770  0.56327139
       [7]  0.14398320  0.58950458 -0.77571740 -0.13483416 -0.39627028  0.44146882

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                           
       [1,] "sp|P39708|DHE5_YEAST" "AANLGGVAVSGLEMAQNSQK"         
       [2,] "sp|P07262|DHE4_YEAST" "AANLGGVAVSGLEMAQNSQR"         
       [3,] "sp|P07262|DHE4_YEAST" "ALVAQGVK"                     
       [4,] "sp|P39708|DHE5_YEAST" "FHPSVNLSILK"                  
       [5,] "sp|P07262|DHE4_YEAST" "FHPSVNLSILK"                  
       [6,] "sp|P07262|DHE4_YEAST" "FIAEGSNMGSTPEAIAVFETAR"       
       [7,] "sp|P39708|DHE5_YEAST" "FLGFEQIFK"                    
       [8,] "sp|P07262|DHE4_YEAST" "FLGFEQIFK"                    
       [9,] "sp|P07262|DHE4_YEAST" "GANIASFIK"                    
      [10,] "sp|P39708|DHE5_YEAST" "GANIASFVMVADAMLDQGDVF"        
      [11,] "sp|P39708|DHE5_YEAST" "GCIISETGITSEQIHDIASAK"        
      [12,] "sp|P07262|DHE4_YEAST" "GCIISETGITSEQVADISSAK"        
      [13,] "sp|P39708|DHE5_YEAST" "GGLCVDLK"                     
      [14,] "sp|P07262|DHE4_YEAST" "GGLCVDLK"                     
      [15,] "sp|P07262|DHE4_YEAST" "HIGQDTDVPAGDIGVGGR"           
      [16,] "sp|P39708|DHE5_YEAST" "ICYAFMR"                      
      [17,] "sp|P07262|DHE4_YEAST" "ICYAFMR"                      
      [18,] "sp|P07262|DHE4_YEAST" "IMINCFNECIDYAK"               
      [19,] "sp|P07262|DHE4_YEAST" "ITWTSER"                      
      [20,] "sp|P39708|DHE5_YEAST" "NSWEGVLTGK"                   
      [21,] "sp|P07262|DHE4_YEAST" "NSWEGVLTGK"                   
      [22,] "sp|P07262|DHE4_YEAST" "SEPEFQQAYEEVVSSLEDSTLFEQHPEYR"
      [23,] "sp|P39708|DHE5_YEAST" "SLEEIVDEYSTFSESK"             
      [24,] "sp|P07262|DHE4_YEAST" "SLEQIVNEYSTFSENK"             
      [25,] "sp|P07262|DHE4_YEAST" "STATGPSEAVWYGPPK"             
      [26,] "sp|P07262|DHE4_YEAST" "VDIALPCATQNEVSGEEAK"          
      [27,] "sp|P07262|DHE4_YEAST" "VIELGGTVVSLSDSK"              
      [28,] "sp|P39708|DHE5_YEAST" "VLPIVSVPER"                   
      [29,] "sp|P07262|DHE4_YEAST" "VLPIVSVPER"                   
      [30,] "sp|P07262|DHE4_YEAST" "VQYIAGARPWTHVQK"              
      [31,] "sp|P39708|DHE5_YEAST" "VTISGSGNVAQYAALK"             
      [32,] "sp|P07262|DHE4_YEAST" "VTISGSGNVAQYAALK"             
      [33,] "sp|P07262|DHE4_YEAST" "VTWENDKGEQEVAQGYR"            
      [34,] "sp|P39708|DHE5_YEAST" "VTWENDNGEQEVAQGYR"            

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA  0.34485174 -0.18344358  0.22632477  0.09075235
       [7]  0.19377265  0.19368206  0.14552313  0.50823246 -0.02677579 -0.38078692
      [13]  0.12241773  0.10749690  0.37479644 -0.07557933 -0.62041663  0.28544390
      [19]  0.16716561  0.03865488 -0.35837351 -0.03775721  0.02416326 -0.14626485
      [25] -0.18176647 -0.56307667 -0.07056382  0.26938130 -0.46539090

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                          
       [1,] "sp|P00330|ADH1_YEAST" "ANELLINVK"                   
       [2,] "sp|P00331|ADH2_YEAST" "ANGTVVLVGLPAGAK"             
       [3,] "sp|P00330|ADH1_YEAST" "ATDGGAHGVINVSVSEAAIEASTR"    
       [4,] "sp|P00330|ADH1_YEAST" "CCSDVFNQVVK"                 
       [5,] "sp|P07246|ADH3_YEAST" "DIPVPEPKPNEILINVK"           
       [6,] "sp|P00330|ADH1_YEAST" "DIVGAVLK"                    
       [7,] "sp|P00330|ADH1_YEAST" "EALDFFAR"                    
       [8,] "sp|P00331|ADH2_YEAST" "EALDFFAR"                    
       [9,] "sp|P38113|ADH5_YEAST" "EALDFFAR"                    
      [10,] "sp|P07246|ADH3_YEAST" "EALDFFSR"                    
      [11,] "sp|P00330|ADH1_YEAST" "EKDIVGAVLK"                  
      [12,] "sp|P07246|ADH3_YEAST" "GVIFYENK"                    
      [13,] "sp|P00330|ADH1_YEAST" "GVIFYESHGK"                  
      [14,] "sp|P00330|ADH1_YEAST" "IGDYAGIK"                    
      [15,] "sp|P00331|ADH2_YEAST" "IGDYAGIK"                    
      [16,] "sp|P07246|ADH3_YEAST" "IQQGTDLAEVAPILCAGVTVYK"      
      [17,] "sp|P07246|ADH3_YEAST" "IVGLSELPK"                   
      [18,] "sp|P00330|ADH1_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [19,] "sp|P00331|ADH2_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [20,] "sp|P07246|ADH3_YEAST" "NMVSDIQEATK"                 
      [21,] "sp|P00330|ADH1_YEAST" "SANLMAGHWVAISGAAGGLGSLAVQYAK"
      [22,] "sp|P00330|ADH1_YEAST" "SIGGEVFIDFTK"                
      [23,] "sp|P00330|ADH1_YEAST" "SIPETQK"                     
      [24,] "sp|P00331|ADH2_YEAST" "SIPETQK"                     
      [25,] "sp|P00330|ADH1_YEAST" "SISIVGSYVGNR"                
      [26,] "sp|P00331|ADH2_YEAST" "SISIVGSYVGNR"                
      [27,] "sp|P00330|ADH1_YEAST" "VLGIDGGEGK"                  
      [28,] "sp|P00330|ADH1_YEAST" "VLGIDGGEGKEELFR"             
      [29,] "sp|P00331|ADH2_YEAST" "VVGLSSLPEIYEK"               
      [30,] "sp|P00330|ADH1_YEAST" "VVGLSTLPEIYEK"               
      [31,] "sp|P00330|ADH1_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       
      [32,] "sp|P07246|ADH3_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA          NA          NA -0.09213118  0.13908309
       [7] -0.63089964  0.18150780  0.10852732  0.01636308  0.02228366  0.07937947
      [13] -0.26744811  0.12388694  0.21514238 -0.06718061  0.10167238  0.20413832
      [19] -0.21856356 -0.60813279  0.02883078  0.13946039  0.54554499 -0.06322847
      [25] -0.94480733  0.22201082 -0.10181185 -0.03660714  0.02723379

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 14

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA  0.942555857  0.191337293  0.473638399  0.281388730
       [6]  0.004830101 -0.143762647 -0.037945653 -0.742235799  0.143241202
      [11] -0.272304407 -0.294216237 -0.070510314 -0.120199738 -0.304583227

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 10

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA  0.01301067 -0.39541482  0.25095770  0.56327139
       [7]  0.14398320  0.58950458 -0.77571740 -0.13483416 -0.39627028  0.44146882

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 27

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA  0.34485174 -0.18344358  0.22632477  0.09075235
       [7]  0.19377265  0.19368206  0.14552313  0.50823246 -0.02677579 -0.38078692
      [13]  0.12241773  0.10749690  0.37479644 -0.07557933 -0.62041663  0.28544390
      [19]  0.16716561  0.03865488 -0.35837351 -0.03775721  0.02416326 -0.14626485
      [25] -0.18176647 -0.56307667 -0.07056382  0.26938130 -0.46539090

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 25

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA          NA          NA -0.09213118  0.13908309
       [7] -0.63089964  0.18150780  0.10852732  0.01636308  0.02228366  0.07937947
      [13] -0.26744811  0.12388694  0.21514238 -0.06718061  0.10167238  0.20413832
      [19] -0.21856356 -0.60813279  0.02883078  0.13946039  0.54554499 -0.06322847
      [25] -0.94480733  0.22201082 -0.10181185 -0.03660714  0.02723379

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                          
       [1,] "sp|P09938|RIR2_YEAST" "AAADALSDLEIKDSK"             
       [2,] "sp|P09938|RIR2_YEAST" "AEWALR"                      
       [3,] "sp|P09938|RIR2_YEAST" "DEGLHTDFACLLFAHLK"           
       [4,] "sp|P09938|RIR2_YEAST" "DIHDWNNR"                    
       [5,] "sp|P09938|RIR2_YEAST" "ESEFLFNAIHTIPEIGEK"          
       [6,] "sp|P09938|RIR2_YEAST" "GMMPGLTFSNELICR"             
       [7,] "sp|P09938|RIR2_YEAST" "IVTEAVEIEQR"                 
       [8,] "sp|P09938|RIR2_YEAST" "LLVAFGNK"                    
       [9,] "sp|P09938|RIR2_YEAST" "NKPDPAIVEK"                  
      [10,] "sp|P09938|RIR2_YEAST" "TNFFEK"                      
      [11,] "sp|P09938|RIR2_YEAST" "TVLFPIK"                     
      [12,] "sp|P09938|RIR2_YEAST" "VENPFDFMENISLAGK"            
      [13,] "sp|P09938|RIR2_YEAST" "WIQDADALFGER"                
      [14,] "sp|P09938|RIR2_YEAST" "YFLDALPVALLGMNADLMNQYVEFVADR"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA  0.111939658  0.114224717  0.834212169  0.166191535
       [6] -0.079398236 -0.063350589  0.126814801  0.007996683  0.008965827
      [11]  0.054036142 -0.041525488 -0.164347183 -0.036941729 -0.165674779

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                    [,2]                           
       [1,] "sp|P40212|RL13B_YEAST" "AAGLTAAYAR"                   
       [2,] "sp|Q12690|RL13A_YEAST" "AAGLTAAYAR"                   
       [3,] "sp|P40212|RL13B_YEAST" "APEAEQVLSAAATFPIAQPATDVEAR"   
       [4,] "sp|Q12690|RL13A_YEAST" "APEAEQVLSAAATFPIAQPATDVEAR"   
       [5,] "sp|P40212|RL13B_YEAST" "AVQDNGESAFR"                  
       [6,] "sp|Q12690|RL13A_YEAST" "AVQDNGESAFR"                  
       [7,] "sp|P40212|RL13B_YEAST" "DGKAPEAEQVLSAAATFPIAQPATDVEAR"
       [8,] "sp|P40212|RL13B_YEAST" "GFTLAEVK"                     
       [9,] "sp|Q12690|RL13A_YEAST" "GFTLAEVK"                     
      [10,] "sp|P40212|RL13B_YEAST" "IAPRPLDLLRPVVR"               
      [11,] "sp|Q12690|RL13A_YEAST" "IAPRPLDLLRPVVR"               
      [12,] "sp|P40212|RL13B_YEAST" "IIVFPR"                       
      [13,] "sp|Q12690|RL13A_YEAST" "IIVFPR"                       
      [14,] "sp|P40212|RL13B_YEAST" "NQEIFDANVQR"                  
      [15,] "sp|Q12690|RL13A_YEAST" "NQEIFDANVQR"                  
      [16,] "sp|P40212|RL13B_YEAST" "TIGIAVDHR"                    
      [17,] "sp|Q12690|RL13A_YEAST" "TIGIAVDHR"                    
      [18,] "sp|P40212|RL13B_YEAST" "VHFDQAGK"                     
      [19,] "sp|Q12690|RL13A_YEAST" "VHFDQAGK"                     

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA  0.147734012  0.088359391  0.098808301
       [6] -0.357043272 -0.012065202 -0.111163879  0.149823716  0.001685123
      [11]  0.060299082  0.148383264

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                           
       [1,] "sp|P39708|DHE5_YEAST" "AANLGGVAVSGLEMAQNSQK"         
       [2,] "sp|P07262|DHE4_YEAST" "AANLGGVAVSGLEMAQNSQR"         
       [3,] "sp|P07262|DHE4_YEAST" "ALVAQGVK"                     
       [4,] "sp|P39708|DHE5_YEAST" "FHPSVNLSILK"                  
       [5,] "sp|P07262|DHE4_YEAST" "FHPSVNLSILK"                  
       [6,] "sp|P07262|DHE4_YEAST" "FIAEGSNMGSTPEAIAVFETAR"       
       [7,] "sp|P39708|DHE5_YEAST" "FLGFEQIFK"                    
       [8,] "sp|P07262|DHE4_YEAST" "FLGFEQIFK"                    
       [9,] "sp|P07262|DHE4_YEAST" "GANIASFIK"                    
      [10,] "sp|P39708|DHE5_YEAST" "GANIASFVMVADAMLDQGDVF"        
      [11,] "sp|P39708|DHE5_YEAST" "GCIISETGITSEQIHDIASAK"        
      [12,] "sp|P07262|DHE4_YEAST" "GCIISETGITSEQVADISSAK"        
      [13,] "sp|P39708|DHE5_YEAST" "GGLCVDLK"                     
      [14,] "sp|P07262|DHE4_YEAST" "GGLCVDLK"                     
      [15,] "sp|P07262|DHE4_YEAST" "HIGQDTDVPAGDIGVGGR"           
      [16,] "sp|P39708|DHE5_YEAST" "ICYAFMR"                      
      [17,] "sp|P07262|DHE4_YEAST" "ICYAFMR"                      
      [18,] "sp|P07262|DHE4_YEAST" "IMINCFNECIDYAK"               
      [19,] "sp|P07262|DHE4_YEAST" "ITWTSER"                      
      [20,] "sp|P39708|DHE5_YEAST" "NSWEGVLTGK"                   
      [21,] "sp|P07262|DHE4_YEAST" "NSWEGVLTGK"                   
      [22,] "sp|P07262|DHE4_YEAST" "SEPEFQQAYEEVVSSLEDSTLFEQHPEYR"
      [23,] "sp|P39708|DHE5_YEAST" "SLEEIVDEYSTFSESK"             
      [24,] "sp|P07262|DHE4_YEAST" "SLEQIVNEYSTFSENK"             
      [25,] "sp|P07262|DHE4_YEAST" "STATGPSEAVWYGPPK"             
      [26,] "sp|P07262|DHE4_YEAST" "VDIALPCATQNEVSGEEAK"          
      [27,] "sp|P07262|DHE4_YEAST" "VIELGGTVVSLSDSK"              
      [28,] "sp|P39708|DHE5_YEAST" "VLPIVSVPER"                   
      [29,] "sp|P07262|DHE4_YEAST" "VLPIVSVPER"                   
      [30,] "sp|P07262|DHE4_YEAST" "VQYIAGARPWTHVQK"              
      [31,] "sp|P39708|DHE5_YEAST" "VTISGSGNVAQYAALK"             
      [32,] "sp|P07262|DHE4_YEAST" "VTISGSGNVAQYAALK"             
      [33,] "sp|P07262|DHE4_YEAST" "VTWENDKGEQEVAQGYR"            
      [34,] "sp|P39708|DHE5_YEAST" "VTWENDNGEQEVAQGYR"            

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA  0.240123495  0.024851220  0.318890663
       [6]  0.521461082  0.150682350 -0.122716431  0.015260966 -0.244874686
      [11] -0.151223269 -0.104874324  0.105964538  0.092190596 -0.277878261
      [16] -0.259936809  0.116524754  0.103673570  0.410506928 -0.109239850
      [21] -0.118725941  0.079106057  0.119379251 -0.021467326  0.198809832
      [26] -0.352768572 -0.008653108 -0.103829079 -0.534700083

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                          
       [1,] "sp|P00330|ADH1_YEAST" "ANELLINVK"                   
       [2,] "sp|P00331|ADH2_YEAST" "ANGTVVLVGLPAGAK"             
       [3,] "sp|P00330|ADH1_YEAST" "ATDGGAHGVINVSVSEAAIEASTR"    
       [4,] "sp|P00330|ADH1_YEAST" "CCSDVFNQVVK"                 
       [5,] "sp|P07246|ADH3_YEAST" "DIPVPEPKPNEILINVK"           
       [6,] "sp|P00330|ADH1_YEAST" "DIVGAVLK"                    
       [7,] "sp|P00330|ADH1_YEAST" "EALDFFAR"                    
       [8,] "sp|P00331|ADH2_YEAST" "EALDFFAR"                    
       [9,] "sp|P38113|ADH5_YEAST" "EALDFFAR"                    
      [10,] "sp|P07246|ADH3_YEAST" "EALDFFSR"                    
      [11,] "sp|P00330|ADH1_YEAST" "EKDIVGAVLK"                  
      [12,] "sp|P07246|ADH3_YEAST" "GVIFYENK"                    
      [13,] "sp|P00330|ADH1_YEAST" "GVIFYESHGK"                  
      [14,] "sp|P00330|ADH1_YEAST" "IGDYAGIK"                    
      [15,] "sp|P00331|ADH2_YEAST" "IGDYAGIK"                    
      [16,] "sp|P07246|ADH3_YEAST" "IQQGTDLAEVAPILCAGVTVYK"      
      [17,] "sp|P07246|ADH3_YEAST" "IVGLSELPK"                   
      [18,] "sp|P00330|ADH1_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [19,] "sp|P00331|ADH2_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [20,] "sp|P07246|ADH3_YEAST" "NMVSDIQEATK"                 
      [21,] "sp|P00330|ADH1_YEAST" "SANLMAGHWVAISGAAGGLGSLAVQYAK"
      [22,] "sp|P00330|ADH1_YEAST" "SIGGEVFIDFTK"                
      [23,] "sp|P00330|ADH1_YEAST" "SIPETQK"                     
      [24,] "sp|P00331|ADH2_YEAST" "SIPETQK"                     
      [25,] "sp|P00330|ADH1_YEAST" "SISIVGSYVGNR"                
      [26,] "sp|P00331|ADH2_YEAST" "SISIVGSYVGNR"                
      [27,] "sp|P00330|ADH1_YEAST" "VLGIDGGEGK"                  
      [28,] "sp|P00330|ADH1_YEAST" "VLGIDGGEGKEELFR"             
      [29,] "sp|P00331|ADH2_YEAST" "VVGLSSLPEIYEK"               
      [30,] "sp|P00330|ADH1_YEAST" "VVGLSTLPEIYEK"               
      [31,] "sp|P00330|ADH1_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       
      [32,] "sp|P07246|ADH3_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA           NA           NA  0.070723529
       [6]  0.079567338 -0.271105299  0.003768771  0.238051007  0.029163650
      [11]  0.040454258 -0.081834924 -0.194013837 -0.047158244 -0.007118395
      [16]  0.243742967 -0.048948307  0.021680998 -0.204002800  0.060607938
      [21]  0.350490567  0.056808858 -0.487930265  0.052001320  0.285153449
      [26] -0.113191177  0.104953674 -0.026544879  0.057434838

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 14

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA  0.111939658  0.114224717  0.834212169  0.166191535
       [6] -0.079398236 -0.063350589  0.126814801  0.007996683  0.008965827
      [11]  0.054036142 -0.041525488 -0.164347183 -0.036941729 -0.165674779

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 10

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA  0.147734012  0.088359391  0.098808301
       [6] -0.357043272 -0.012065202 -0.111163879  0.149823716  0.001685123
      [11]  0.060299082  0.148383264

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 27

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA  0.240123495  0.024851220  0.318890663
       [6]  0.521461082  0.150682350 -0.122716431  0.015260966 -0.244874686
      [11] -0.151223269 -0.104874324  0.105964538  0.092190596 -0.277878261
      [16] -0.259936809  0.116524754  0.103673570  0.410506928 -0.109239850
      [21] -0.118725941  0.079106057  0.119379251 -0.021467326  0.198809832
      [26] -0.352768572 -0.008653108 -0.103829079 -0.534700083

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 25

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA           NA           NA  0.070723529
       [6]  0.079567338 -0.271105299  0.003768771  0.238051007  0.029163650
      [11]  0.040454258 -0.081834924 -0.194013837 -0.047158244 -0.007118395
      [16]  0.243742967 -0.048948307  0.021680998 -0.204002800  0.060607938
      [21]  0.350490567  0.056808858 -0.487930265  0.052001320  0.285153449
      [26] -0.113191177  0.104953674 -0.026544879  0.057434838

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                          
       [1,] "sp|P09938|RIR2_YEAST" "AAADALSDLEIKDSK"             
       [2,] "sp|P09938|RIR2_YEAST" "AEWALR"                      
       [3,] "sp|P09938|RIR2_YEAST" "DEGLHTDFACLLFAHLK"           
       [4,] "sp|P09938|RIR2_YEAST" "DIHDWNNR"                    
       [5,] "sp|P09938|RIR2_YEAST" "ESEFLFNAIHTIPEIGEK"          
       [6,] "sp|P09938|RIR2_YEAST" "GMMPGLTFSNELICR"             
       [7,] "sp|P09938|RIR2_YEAST" "IVTEAVEIEQR"                 
       [8,] "sp|P09938|RIR2_YEAST" "LLVAFGNK"                    
       [9,] "sp|P09938|RIR2_YEAST" "NKPDPAIVEK"                  
      [10,] "sp|P09938|RIR2_YEAST" "TNFFEK"                      
      [11,] "sp|P09938|RIR2_YEAST" "TVLFPIK"                     
      [12,] "sp|P09938|RIR2_YEAST" "VENPFDFMENISLAGK"            
      [13,] "sp|P09938|RIR2_YEAST" "WIQDADALFGER"                
      [14,] "sp|P09938|RIR2_YEAST" "YFLDALPVALLGMNADLMNQYVEFVADR"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA  0.067985566  0.008548272  0.714583859 -1.756814764
       [6] -0.074776715  0.149209984 -0.025286382 -0.237088193 -0.042675710
      [11] -0.061904228 -0.235351754 -0.115624288  0.172954743 -0.180940859

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                    [,2]                           
       [1,] "sp|P40212|RL13B_YEAST" "AAGLTAAYAR"                   
       [2,] "sp|Q12690|RL13A_YEAST" "AAGLTAAYAR"                   
       [3,] "sp|P40212|RL13B_YEAST" "APEAEQVLSAAATFPIAQPATDVEAR"   
       [4,] "sp|Q12690|RL13A_YEAST" "APEAEQVLSAAATFPIAQPATDVEAR"   
       [5,] "sp|P40212|RL13B_YEAST" "AVQDNGESAFR"                  
       [6,] "sp|Q12690|RL13A_YEAST" "AVQDNGESAFR"                  
       [7,] "sp|P40212|RL13B_YEAST" "DGKAPEAEQVLSAAATFPIAQPATDVEAR"
       [8,] "sp|P40212|RL13B_YEAST" "GFTLAEVK"                     
       [9,] "sp|Q12690|RL13A_YEAST" "GFTLAEVK"                     
      [10,] "sp|P40212|RL13B_YEAST" "IAPRPLDLLRPVVR"               
      [11,] "sp|Q12690|RL13A_YEAST" "IAPRPLDLLRPVVR"               
      [12,] "sp|P40212|RL13B_YEAST" "IIVFPR"                       
      [13,] "sp|Q12690|RL13A_YEAST" "IIVFPR"                       
      [14,] "sp|P40212|RL13B_YEAST" "NQEIFDANVQR"                  
      [15,] "sp|Q12690|RL13A_YEAST" "NQEIFDANVQR"                  
      [16,] "sp|P40212|RL13B_YEAST" "TIGIAVDHR"                    
      [17,] "sp|Q12690|RL13A_YEAST" "TIGIAVDHR"                    
      [18,] "sp|P40212|RL13B_YEAST" "VHFDQAGK"                     
      [19,] "sp|Q12690|RL13A_YEAST" "VHFDQAGK"                     

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA -0.06266107 -0.00938573  0.14562423 -0.19616627
       [7]  0.10162963 -0.41162886  0.19918714  0.04051594 -0.08897282 -0.01023307

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                           
       [1,] "sp|P39708|DHE5_YEAST" "AANLGGVAVSGLEMAQNSQK"         
       [2,] "sp|P07262|DHE4_YEAST" "AANLGGVAVSGLEMAQNSQR"         
       [3,] "sp|P07262|DHE4_YEAST" "ALVAQGVK"                     
       [4,] "sp|P39708|DHE5_YEAST" "FHPSVNLSILK"                  
       [5,] "sp|P07262|DHE4_YEAST" "FHPSVNLSILK"                  
       [6,] "sp|P07262|DHE4_YEAST" "FIAEGSNMGSTPEAIAVFETAR"       
       [7,] "sp|P39708|DHE5_YEAST" "FLGFEQIFK"                    
       [8,] "sp|P07262|DHE4_YEAST" "FLGFEQIFK"                    
       [9,] "sp|P07262|DHE4_YEAST" "GANIASFIK"                    
      [10,] "sp|P39708|DHE5_YEAST" "GANIASFVMVADAMLDQGDVF"        
      [11,] "sp|P39708|DHE5_YEAST" "GCIISETGITSEQIHDIASAK"        
      [12,] "sp|P07262|DHE4_YEAST" "GCIISETGITSEQVADISSAK"        
      [13,] "sp|P39708|DHE5_YEAST" "GGLCVDLK"                     
      [14,] "sp|P07262|DHE4_YEAST" "GGLCVDLK"                     
      [15,] "sp|P07262|DHE4_YEAST" "HIGQDTDVPAGDIGVGGR"           
      [16,] "sp|P39708|DHE5_YEAST" "ICYAFMR"                      
      [17,] "sp|P07262|DHE4_YEAST" "ICYAFMR"                      
      [18,] "sp|P07262|DHE4_YEAST" "IMINCFNECIDYAK"               
      [19,] "sp|P07262|DHE4_YEAST" "ITWTSER"                      
      [20,] "sp|P39708|DHE5_YEAST" "NSWEGVLTGK"                   
      [21,] "sp|P07262|DHE4_YEAST" "NSWEGVLTGK"                   
      [22,] "sp|P07262|DHE4_YEAST" "SEPEFQQAYEEVVSSLEDSTLFEQHPEYR"
      [23,] "sp|P39708|DHE5_YEAST" "SLEEIVDEYSTFSESK"             
      [24,] "sp|P07262|DHE4_YEAST" "SLEQIVNEYSTFSENK"             
      [25,] "sp|P07262|DHE4_YEAST" "STATGPSEAVWYGPPK"             
      [26,] "sp|P07262|DHE4_YEAST" "VDIALPCATQNEVSGEEAK"          
      [27,] "sp|P07262|DHE4_YEAST" "VIELGGTVVSLSDSK"              
      [28,] "sp|P39708|DHE5_YEAST" "VLPIVSVPER"                   
      [29,] "sp|P07262|DHE4_YEAST" "VLPIVSVPER"                   
      [30,] "sp|P07262|DHE4_YEAST" "VQYIAGARPWTHVQK"              
      [31,] "sp|P39708|DHE5_YEAST" "VTISGSGNVAQYAALK"             
      [32,] "sp|P07262|DHE4_YEAST" "VTISGSGNVAQYAALK"             
      [33,] "sp|P07262|DHE4_YEAST" "VTWENDKGEQEVAQGYR"            
      [34,] "sp|P39708|DHE5_YEAST" "VTWENDNGEQEVAQGYR"            

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA  0.360632543  0.191592940  0.069264628
       [6]  0.357116431  0.047761246  0.022189996  0.104338487 -0.018934983
      [11] -0.158458914 -0.078285974  0.013855982  0.034193021 -0.085904775
      [16]  0.014579584 -0.070022059  0.003081272  0.514644400 -0.001940633
      [21] -0.109404832  0.182301634  0.215443910 -0.060726809  0.116976595
      [26]  0.356825909 -0.239825860 -0.040603877 -0.139176241

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                          
       [1,] "sp|P00330|ADH1_YEAST" "ANELLINVK"                   
       [2,] "sp|P00331|ADH2_YEAST" "ANGTVVLVGLPAGAK"             
       [3,] "sp|P00330|ADH1_YEAST" "ATDGGAHGVINVSVSEAAIEASTR"    
       [4,] "sp|P00330|ADH1_YEAST" "CCSDVFNQVVK"                 
       [5,] "sp|P07246|ADH3_YEAST" "DIPVPEPKPNEILINVK"           
       [6,] "sp|P00330|ADH1_YEAST" "DIVGAVLK"                    
       [7,] "sp|P00330|ADH1_YEAST" "EALDFFAR"                    
       [8,] "sp|P00331|ADH2_YEAST" "EALDFFAR"                    
       [9,] "sp|P38113|ADH5_YEAST" "EALDFFAR"                    
      [10,] "sp|P07246|ADH3_YEAST" "EALDFFSR"                    
      [11,] "sp|P00330|ADH1_YEAST" "EKDIVGAVLK"                  
      [12,] "sp|P07246|ADH3_YEAST" "GVIFYENK"                    
      [13,] "sp|P00330|ADH1_YEAST" "GVIFYESHGK"                  
      [14,] "sp|P00330|ADH1_YEAST" "IGDYAGIK"                    
      [15,] "sp|P00331|ADH2_YEAST" "IGDYAGIK"                    
      [16,] "sp|P07246|ADH3_YEAST" "IQQGTDLAEVAPILCAGVTVYK"      
      [17,] "sp|P07246|ADH3_YEAST" "IVGLSELPK"                   
      [18,] "sp|P00330|ADH1_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [19,] "sp|P00331|ADH2_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [20,] "sp|P07246|ADH3_YEAST" "NMVSDIQEATK"                 
      [21,] "sp|P00330|ADH1_YEAST" "SANLMAGHWVAISGAAGGLGSLAVQYAK"
      [22,] "sp|P00330|ADH1_YEAST" "SIGGEVFIDFTK"                
      [23,] "sp|P00330|ADH1_YEAST" "SIPETQK"                     
      [24,] "sp|P00331|ADH2_YEAST" "SIPETQK"                     
      [25,] "sp|P00330|ADH1_YEAST" "SISIVGSYVGNR"                
      [26,] "sp|P00331|ADH2_YEAST" "SISIVGSYVGNR"                
      [27,] "sp|P00330|ADH1_YEAST" "VLGIDGGEGK"                  
      [28,] "sp|P00330|ADH1_YEAST" "VLGIDGGEGKEELFR"             
      [29,] "sp|P00331|ADH2_YEAST" "VVGLSSLPEIYEK"               
      [30,] "sp|P00330|ADH1_YEAST" "VVGLSTLPEIYEK"               
      [31,] "sp|P00330|ADH1_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       
      [32,] "sp|P07246|ADH3_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA           NA           NA -0.165942403
       [6] -0.009944959 -0.402737246  0.115159014  0.304261246 -0.090984073
      [11]  0.082660029  0.187851888 -0.261168113 -0.016081309 -0.088188097
      [16]  0.019101663 -0.006168545  0.219456480 -0.045053522  0.069013948
      [21]  0.118699361  0.037229437  0.068071219 -0.013481381 -0.043138486
      [26]  0.028660015 -0.025245880  0.009888782  0.606720613

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 14

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA  0.067985566  0.008548272  0.714583859 -1.756814764
       [6] -0.074776715  0.149209984 -0.025286382 -0.237088193 -0.042675710
      [11] -0.061904228 -0.235351754 -0.115624288  0.172954743 -0.180940859

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 10

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA -0.06266107 -0.00938573  0.14562423 -0.19616627
       [7]  0.10162963 -0.41162886  0.19918714  0.04051594 -0.08897282 -0.01023307

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 27

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA  0.360632543  0.191592940  0.069264628
       [6]  0.357116431  0.047761246  0.022189996  0.104338487 -0.018934983
      [11] -0.158458914 -0.078285974  0.013855982  0.034193021 -0.085904775
      [16]  0.014579584 -0.070022059  0.003081272  0.514644400 -0.001940633
      [21] -0.109404832  0.182301634  0.215443910 -0.060726809  0.116976595
      [26]  0.356825909 -0.239825860 -0.040603877 -0.139176241

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 25

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA           NA           NA -0.165942403
       [6] -0.009944959 -0.402737246  0.115159014  0.304261246 -0.090984073
      [11]  0.082660029  0.187851888 -0.261168113 -0.016081309 -0.088188097
      [16]  0.019101663 -0.006168545  0.219456480 -0.045053522  0.069013948
      [21]  0.118699361  0.037229437  0.068071219 -0.013481381 -0.043138486
      [26]  0.028660015 -0.025245880  0.009888782  0.606720613

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                          
       [1,] "sp|P09938|RIR2_YEAST" "AAADALSDLEIKDSK"             
       [2,] "sp|P09938|RIR2_YEAST" "AEWALR"                      
       [3,] "sp|P09938|RIR2_YEAST" "DIHDWNNR"                    
       [4,] "sp|P09938|RIR2_YEAST" "ESEFLFNAIHTIPEIGEK"          
       [5,] "sp|P09938|RIR2_YEAST" "GMMPGLTFSNELICR"             
       [6,] "sp|P09938|RIR2_YEAST" "IVTEAVEIEQR"                 
       [7,] "sp|P09938|RIR2_YEAST" "LLVAFGNK"                    
       [8,] "sp|P09938|RIR2_YEAST" "NKPDPAIVEK"                  
       [9,] "sp|P09938|RIR2_YEAST" "TNFFEK"                      
      [10,] "sp|P09938|RIR2_YEAST" "TVLFPIK"                     
      [11,] "sp|P09938|RIR2_YEAST" "VENPFDFMENISLAGK"            
      [12,] "sp|P09938|RIR2_YEAST" "WIQDADALFGER"                
      [13,] "sp|P09938|RIR2_YEAST" "YFLDALPVALLGMNADLMNQYVEFVADR"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA  0.840405327  0.116390221  0.034122563 -0.209968681
       [6] -0.072845446 -0.003357929 -0.782032752  0.197397630 -0.110266104
      [11] -0.337944782  0.045976768 -0.173114701 -0.212709022

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                    [,2]                           
       [1,] "sp|P40212|RL13B_YEAST" "AAGLTAAYAR"                   
       [2,] "sp|Q12690|RL13A_YEAST" "AAGLTAAYAR"                   
       [3,] "sp|P40212|RL13B_YEAST" "APEAEQVLSAAATFPIAQPATDVEAR"   
       [4,] "sp|Q12690|RL13A_YEAST" "APEAEQVLSAAATFPIAQPATDVEAR"   
       [5,] "sp|P40212|RL13B_YEAST" "AVQDNGESAFR"                  
       [6,] "sp|Q12690|RL13A_YEAST" "AVQDNGESAFR"                  
       [7,] "sp|P40212|RL13B_YEAST" "DGKAPEAEQVLSAAATFPIAQPATDVEAR"
       [8,] "sp|P40212|RL13B_YEAST" "GFTLAEVK"                     
       [9,] "sp|Q12690|RL13A_YEAST" "GFTLAEVK"                     
      [10,] "sp|P40212|RL13B_YEAST" "IAPRPLDLLRPVVR"               
      [11,] "sp|Q12690|RL13A_YEAST" "IAPRPLDLLRPVVR"               
      [12,] "sp|P40212|RL13B_YEAST" "IIVFPR"                       
      [13,] "sp|Q12690|RL13A_YEAST" "IIVFPR"                       
      [14,] "sp|P40212|RL13B_YEAST" "NQEIFDANVQR"                  
      [15,] "sp|Q12690|RL13A_YEAST" "NQEIFDANVQR"                  
      [16,] "sp|P40212|RL13B_YEAST" "TIGIAVDHR"                    
      [17,] "sp|Q12690|RL13A_YEAST" "TIGIAVDHR"                    
      [18,] "sp|P40212|RL13B_YEAST" "VHFDQAGK"                     
      [19,] "sp|Q12690|RL13A_YEAST" "VHFDQAGK"                     

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA  0.03479784 -0.63422266 -0.44659286  0.01836650
       [7]  0.15001875  0.26970881 -0.13678309 -0.01924402 -0.16607214  0.10623702

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                           
       [1,] "sp|P39708|DHE5_YEAST" "AANLGGVAVSGLEMAQNSQK"         
       [2,] "sp|P07262|DHE4_YEAST" "AANLGGVAVSGLEMAQNSQR"         
       [3,] "sp|P07262|DHE4_YEAST" "ALVAQGVK"                     
       [4,] "sp|P39708|DHE5_YEAST" "FHPSVNLSILK"                  
       [5,] "sp|P07262|DHE4_YEAST" "FHPSVNLSILK"                  
       [6,] "sp|P07262|DHE4_YEAST" "FIAEGSNMGSTPEAIAVFETAR"       
       [7,] "sp|P39708|DHE5_YEAST" "FLGFEQIFK"                    
       [8,] "sp|P07262|DHE4_YEAST" "FLGFEQIFK"                    
       [9,] "sp|P07262|DHE4_YEAST" "GANIASFIK"                    
      [10,] "sp|P39708|DHE5_YEAST" "GANIASFVMVADAMLDQGDVF"        
      [11,] "sp|P39708|DHE5_YEAST" "GCIISETGITSEQIHDIASAK"        
      [12,] "sp|P07262|DHE4_YEAST" "GCIISETGITSEQVADISSAK"        
      [13,] "sp|P39708|DHE5_YEAST" "GGLCVDLK"                     
      [14,] "sp|P07262|DHE4_YEAST" "GGLCVDLK"                     
      [15,] "sp|P07262|DHE4_YEAST" "HIGQDTDVPAGDIGVGGR"           
      [16,] "sp|P39708|DHE5_YEAST" "ICYAFMR"                      
      [17,] "sp|P07262|DHE4_YEAST" "ICYAFMR"                      
      [18,] "sp|P07262|DHE4_YEAST" "IMINCFNECIDYAK"               
      [19,] "sp|P07262|DHE4_YEAST" "ITWTSER"                      
      [20,] "sp|P39708|DHE5_YEAST" "NSWEGVLTGK"                   
      [21,] "sp|P07262|DHE4_YEAST" "NSWEGVLTGK"                   
      [22,] "sp|P07262|DHE4_YEAST" "SEPEFQQAYEEVVSSLEDSTLFEQHPEYR"
      [23,] "sp|P39708|DHE5_YEAST" "SLEEIVDEYSTFSESK"             
      [24,] "sp|P07262|DHE4_YEAST" "SLEQIVNEYSTFSENK"             
      [25,] "sp|P07262|DHE4_YEAST" "STATGPSEAVWYGPPK"             
      [26,] "sp|P07262|DHE4_YEAST" "VDIALPCATQNEVSGEEAK"          
      [27,] "sp|P07262|DHE4_YEAST" "VIELGGTVVSLSDSK"              
      [28,] "sp|P39708|DHE5_YEAST" "VLPIVSVPER"                   
      [29,] "sp|P07262|DHE4_YEAST" "VLPIVSVPER"                   
      [30,] "sp|P07262|DHE4_YEAST" "VQYIAGARPWTHVQK"              
      [31,] "sp|P39708|DHE5_YEAST" "VTISGSGNVAQYAALK"             
      [32,] "sp|P07262|DHE4_YEAST" "VTISGSGNVAQYAALK"             
      [33,] "sp|P07262|DHE4_YEAST" "VTWENDKGEQEVAQGYR"            
      [34,] "sp|P39708|DHE5_YEAST" "VTWENDNGEQEVAQGYR"            

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]            NA            NA  0.2534737986  0.3438856627  0.1376682492
       [6]  0.2686097395  0.0535723798  0.1631549138  0.1487584233  0.5048040237
      [11] -0.0206202525  0.0003247414  0.1914329449  0.1058827853  0.2540715227
      [16]  0.0715839437 -0.0243024646  0.1290004143  0.6708864931 -0.1101437044
      [21] -0.0948059632  0.0103997588 -0.0126858870 -0.0850681449 -0.1293403987
      [26]  0.5847865900 -1.0381501336  0.1635293818 -0.2843747810

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                          
       [1,] "sp|P00330|ADH1_YEAST" "ANELLINVK"                   
       [2,] "sp|P00331|ADH2_YEAST" "ANGTVVLVGLPAGAK"             
       [3,] "sp|P00330|ADH1_YEAST" "ATDGGAHGVINVSVSEAAIEASTR"    
       [4,] "sp|P00330|ADH1_YEAST" "CCSDVFNQVVK"                 
       [5,] "sp|P07246|ADH3_YEAST" "DIPVPEPKPNEILINVK"           
       [6,] "sp|P00330|ADH1_YEAST" "DIVGAVLK"                    
       [7,] "sp|P00330|ADH1_YEAST" "EALDFFAR"                    
       [8,] "sp|P00331|ADH2_YEAST" "EALDFFAR"                    
       [9,] "sp|P38113|ADH5_YEAST" "EALDFFAR"                    
      [10,] "sp|P07246|ADH3_YEAST" "EALDFFSR"                    
      [11,] "sp|P00330|ADH1_YEAST" "EKDIVGAVLK"                  
      [12,] "sp|P07246|ADH3_YEAST" "GVIFYENK"                    
      [13,] "sp|P00330|ADH1_YEAST" "GVIFYESHGK"                  
      [14,] "sp|P00330|ADH1_YEAST" "IGDYAGIK"                    
      [15,] "sp|P00331|ADH2_YEAST" "IGDYAGIK"                    
      [16,] "sp|P07246|ADH3_YEAST" "IQQGTDLAEVAPILCAGVTVYK"      
      [17,] "sp|P07246|ADH3_YEAST" "IVGLSELPK"                   
      [18,] "sp|P00330|ADH1_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [19,] "sp|P00331|ADH2_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [20,] "sp|P07246|ADH3_YEAST" "NMVSDIQEATK"                 
      [21,] "sp|P00330|ADH1_YEAST" "SANLMAGHWVAISGAAGGLGSLAVQYAK"
      [22,] "sp|P00330|ADH1_YEAST" "SIGGEVFIDFTK"                
      [23,] "sp|P00330|ADH1_YEAST" "SIPETQK"                     
      [24,] "sp|P00331|ADH2_YEAST" "SIPETQK"                     
      [25,] "sp|P00330|ADH1_YEAST" "SISIVGSYVGNR"                
      [26,] "sp|P00331|ADH2_YEAST" "SISIVGSYVGNR"                
      [27,] "sp|P00330|ADH1_YEAST" "VLGIDGGEGK"                  
      [28,] "sp|P00330|ADH1_YEAST" "VLGIDGGEGKEELFR"             
      [29,] "sp|P00331|ADH2_YEAST" "VVGLSSLPEIYEK"               
      [30,] "sp|P00330|ADH1_YEAST" "VVGLSTLPEIYEK"               
      [31,] "sp|P00330|ADH1_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       
      [32,] "sp|P07246|ADH3_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA           NA           NA -0.061532659
       [6]  0.057541179 -0.610953727  0.166435857 -0.399195511  0.067249138
      [11]  0.044749865  0.150119735 -0.132210732  0.137446952  0.044062637
      [16] -0.031481565  0.210754820  0.178952464 -0.230811399 -0.030650969
      [21]  0.099617676 -0.109113758 -0.071731841 -0.067007223 -0.944393615
      [26]  0.135148664 -0.117533079  0.006010453  0.138091895

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 13

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA  0.840405327  0.116390221  0.034122563 -0.209968681
       [6] -0.072845446 -0.003357929 -0.782032752  0.197397630 -0.110266104
      [11] -0.337944782  0.045976768 -0.173114701 -0.212709022

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 10

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA  0.03479784 -0.63422266 -0.44659286  0.01836650
       [7]  0.15001875  0.26970881 -0.13678309 -0.01924402 -0.16607214  0.10623702

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 27

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]            NA            NA  0.2534737986  0.3438856627  0.1376682492
       [6]  0.2686097395  0.0535723798  0.1631549138  0.1487584233  0.5048040237
      [11] -0.0206202525  0.0003247414  0.1914329449  0.1058827853  0.2540715227
      [16]  0.0715839437 -0.0243024646  0.1290004143  0.6708864931 -0.1101437044
      [21] -0.0948059632  0.0103997588 -0.0126858870 -0.0850681449 -0.1293403987
      [26]  0.5847865900 -1.0381501336  0.1635293818 -0.2843747810

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 25

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA           NA           NA -0.061532659
       [6]  0.057541179 -0.610953727  0.166435857 -0.399195511  0.067249138
      [11]  0.044749865  0.150119735 -0.132210732  0.137446952  0.044062637
      [16] -0.031481565  0.210754820  0.178952464 -0.230811399 -0.030650969
      [21]  0.099617676 -0.109113758 -0.071731841 -0.067007223 -0.944393615
      [26]  0.135148664 -0.117533079  0.006010453  0.138091895

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

