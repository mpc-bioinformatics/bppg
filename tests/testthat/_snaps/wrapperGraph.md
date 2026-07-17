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

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA  0.19668773  2.62855629 -0.22465493  0.61038124
       [7]  0.01046291  0.96942305  0.93449142  0.07345646 -0.60412647

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                           
       [1,] "sp|P39708|DHE5_YEAST" "AANLGGVAVSGLEMAQNSQK"         
       [2,] "sp|P07262|DHE4_YEAST" "AANLGGVAVSGLEMAQNSQR"         
       [3,] "sp|P39708|DHE5_YEAST" "EIGYLFGAYR"                   
       [4,] "sp|P07262|DHE4_YEAST" "EIGYLFGAYR"                   
       [5,] "sp|P39708|DHE5_YEAST" "FHPSVNLSILK"                  
       [6,] "sp|P07262|DHE4_YEAST" "FHPSVNLSILK"                  
       [7,] "sp|P07262|DHE4_YEAST" "FIAEGSNMGSTPEAIAVFETAR"       
       [8,] "sp|P39708|DHE5_YEAST" "FLGFEQIFK"                    
       [9,] "sp|P07262|DHE4_YEAST" "FLGFEQIFK"                    
      [10,] "sp|P07262|DHE4_YEAST" "GANIASFIK"                    
      [11,] "sp|P39708|DHE5_YEAST" "GANIASFVMVADAMLDQGDVF"        
      [12,] "sp|P39708|DHE5_YEAST" "GCIISETGITSEQIHDIASAK"        
      [13,] "sp|P07262|DHE4_YEAST" "GCIISETGITSEQVADISSAK"        
      [14,] "sp|P07262|DHE4_YEAST" "HIGQDTDVPAGDIGVGGR"           
      [15,] "sp|P39708|DHE5_YEAST" "ICYAFMR"                      
      [16,] "sp|P07262|DHE4_YEAST" "ICYAFMR"                      
      [17,] "sp|P07262|DHE4_YEAST" "IMINCFNECIDYAK"               
      [18,] "sp|P07262|DHE4_YEAST" "ITWTSER"                      
      [19,] "sp|P39708|DHE5_YEAST" "NSWEGVLTGK"                   
      [20,] "sp|P07262|DHE4_YEAST" "NSWEGVLTGK"                   
      [21,] "sp|P07262|DHE4_YEAST" "SEPEFQQAYEEVVSSLEDSTLFEQHPEYR"
      [22,] "sp|P39708|DHE5_YEAST" "SLEEIVDEYSTFSESK"             
      [23,] "sp|P07262|DHE4_YEAST" "SLEQIVNEYSTFSENK"             
      [24,] "sp|P07262|DHE4_YEAST" "STATGPSEAVWYGPPK"             
      [25,] "sp|P07262|DHE4_YEAST" "VDIALPCATQNEVSGEEAK"          
      [26,] "sp|P07262|DHE4_YEAST" "VIELGGTVVSLSDSK"              
      [27,] "sp|P39708|DHE5_YEAST" "VLPIVSVPER"                   
      [28,] "sp|P07262|DHE4_YEAST" "VLPIVSVPER"                   
      [29,] "sp|P39708|DHE5_YEAST" "VTISGSGNVAQYAALK"             
      [30,] "sp|P07262|DHE4_YEAST" "VTISGSGNVAQYAALK"             
      [31,] "sp|P07262|DHE4_YEAST" "VTWENDKGEQEVAQGYR"            
      [32,] "sp|P39708|DHE5_YEAST" "VTWENDNGEQEVAQGYR"            

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA  0.09308867 -0.10837982 -0.08378687 -0.69032528
       [7]  0.13659486 -0.19543070  0.11822951 -0.23578507 -0.04787399 -0.62894291
      [13] -0.07947540  0.31576138 -0.14960230 -0.67861068  0.13261953 -0.31428216
      [19]  0.35779330  0.54886191  0.17464817  0.60731745 -0.34997383  0.15436625
      [25]  0.78505065 -0.22981458  0.19616144

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                          
       [1,] "sp|P09938|RIR2_YEAST" "AEWALR"                      
       [2,] "sp|P09938|RIR2_YEAST" "DEGLHTDFACLLFAHLK"           
       [3,] "sp|P09938|RIR2_YEAST" "DIHDWNNR"                    
       [4,] "sp|P09938|RIR2_YEAST" "ELETLREENR"                  
       [5,] "sp|P09938|RIR2_YEAST" "ESEFLFNAIHTIPEIGEK"          
       [6,] "sp|P09938|RIR2_YEAST" "IVTEAVEIEQR"                 
       [7,] "sp|P09938|RIR2_YEAST" "LLVAFGNK"                    
       [8,] "sp|P09938|RIR2_YEAST" "LLVAFGNKK"                   
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
       [1]           NA  0.374014172 -0.766729603  0.657054792 -0.148918279
       [6]  0.188694131 -0.410775258  0.053253590  0.008720694 -0.112580059
      [11]  0.060211627  0.031470669  0.223757471  0.136399210 -0.311319256

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                          
       [1,] "sp|P00330|ADH1_YEAST" "ANELLINVK"                   
       [2,] "sp|P00331|ADH2_YEAST" "ANGTVVLVGLPAGAK"             
       [3,] "sp|P00330|ADH1_YEAST" "ATDGGAHGVINVSVSEAAIEASTR"    
       [4,] "sp|P00330|ADH1_YEAST" "CCSDVFNQVVK"                 
       [5,] "sp|P00331|ADH2_YEAST" "CSSDVFNHVVK"                 
       [6,] "sp|P07246|ADH3_YEAST" "DIPVPEPKPNEILINVK"           
       [7,] "sp|P00330|ADH1_YEAST" "DIVGAVLK"                    
       [8,] "sp|P00330|ADH1_YEAST" "EALDFFAR"                    
       [9,] "sp|P00331|ADH2_YEAST" "EALDFFAR"                    
      [10,] "sp|P38113|ADH5_YEAST" "EALDFFAR"                    
      [11,] "sp|P07246|ADH3_YEAST" "EALDFFSR"                    
      [12,] "sp|P00330|ADH1_YEAST" "EKDIVGAVLK"                  
      [13,] "sp|P07246|ADH3_YEAST" "GVIFYENK"                    
      [14,] "sp|P00330|ADH1_YEAST" "GVIFYESHGK"                  
      [15,] "sp|P00330|ADH1_YEAST" "IGDYAGIK"                    
      [16,] "sp|P00331|ADH2_YEAST" "IGDYAGIK"                    
      [17,] "sp|P07246|ADH3_YEAST" "IQQGTDLAEVAPILCAGVTVYK"      
      [18,] "sp|P07246|ADH3_YEAST" "IVGLSELPK"                   
      [19,] "sp|P00330|ADH1_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [20,] "sp|P00331|ADH2_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [21,] "sp|P07246|ADH3_YEAST" "NMVSDIQEATK"                 
      [22,] "sp|P00330|ADH1_YEAST" "SANLMAGHWVAISGAAGGLGSLAVQYAK"
      [23,] "sp|P00330|ADH1_YEAST" "SIGGEVFIDFTK"                
      [24,] "sp|P00330|ADH1_YEAST" "SIPETQK"                     
      [25,] "sp|P00331|ADH2_YEAST" "SIPETQK"                     
      [26,] "sp|P00330|ADH1_YEAST" "SISIVGSYVGNR"                
      [27,] "sp|P00331|ADH2_YEAST" "SISIVGSYVGNR"                
      [28,] "sp|P00330|ADH1_YEAST" "VLGIDGGEGK"                  
      [29,] "sp|P00330|ADH1_YEAST" "VLGIDGGEGKEELFR"             
      [30,] "sp|P00331|ADH2_YEAST" "VVGLSSLPEIYEK"               
      [31,] "sp|P00330|ADH1_YEAST" "VVGLSTLPEIYEK"               
      [32,] "sp|P00330|ADH1_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       
      [33,] "sp|P07246|ADH3_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA           NA           NA  0.016645753
       [6]  0.237805990 -0.328780267 -0.107871724  0.199970666  0.244710647
      [11]  0.012800875  0.004805219  0.310759179 -0.193894616  0.076038214
      [16] -0.731446490  0.157708335  0.206161441  0.300134275  0.019823451
      [21] -0.561472030 -0.678659245  0.308357062  0.083650408  0.049012247
      [26]  0.142249370 -0.232102345  0.261703095  0.053340916  0.109166744

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 9

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA  0.19668773  2.62855629 -0.22465493  0.61038124
       [7]  0.01046291  0.96942305  0.93449142  0.07345646 -0.60412647

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
       [1]          NA          NA  0.09308867 -0.10837982 -0.08378687 -0.69032528
       [7]  0.13659486 -0.19543070  0.11822951 -0.23578507 -0.04787399 -0.62894291
      [13] -0.07947540  0.31576138 -0.14960230 -0.67861068  0.13261953 -0.31428216
      [19]  0.35779330  0.54886191  0.17464817  0.60731745 -0.34997383  0.15436625
      [25]  0.78505065 -0.22981458  0.19616144

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 14

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA  0.374014172 -0.766729603  0.657054792 -0.148918279
       [6]  0.188694131 -0.410775258  0.053253590  0.008720694 -0.112580059
      [11]  0.060211627  0.031470669  0.223757471  0.136399210 -0.311319256

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 26

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA           NA           NA  0.016645753
       [6]  0.237805990 -0.328780267 -0.107871724  0.199970666  0.244710647
      [11]  0.012800875  0.004805219  0.310759179 -0.193894616  0.076038214
      [16] -0.731446490  0.157708335  0.206161441  0.300134275  0.019823451
      [21] -0.561472030 -0.678659245  0.308357062  0.083650408  0.049012247
      [26]  0.142249370 -0.232102345  0.261703095  0.053340916  0.109166744

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

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
       [1]          NA          NA  0.21847490  2.38974845 -0.92220549  0.06547634
       [7]  0.01649847  0.64962728  1.57342573  0.18904660 -0.37392833 -0.42862844

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                           
       [1,] "sp|P39708|DHE5_YEAST" "AANLGGVAVSGLEMAQNSQK"         
       [2,] "sp|P07262|DHE4_YEAST" "AANLGGVAVSGLEMAQNSQR"         
       [3,] "sp|P39708|DHE5_YEAST" "EIGYLFGAYR"                   
       [4,] "sp|P07262|DHE4_YEAST" "EIGYLFGAYR"                   
       [5,] "sp|P39708|DHE5_YEAST" "FHPSVNLSILK"                  
       [6,] "sp|P07262|DHE4_YEAST" "FHPSVNLSILK"                  
       [7,] "sp|P07262|DHE4_YEAST" "FIAEGSNMGSTPEAIAVFETAR"       
       [8,] "sp|P39708|DHE5_YEAST" "FLGFEQIFK"                    
       [9,] "sp|P07262|DHE4_YEAST" "FLGFEQIFK"                    
      [10,] "sp|P07262|DHE4_YEAST" "GANIASFIK"                    
      [11,] "sp|P39708|DHE5_YEAST" "GANIASFVMVADAMLDQGDVF"        
      [12,] "sp|P39708|DHE5_YEAST" "GCIISETGITSEQIHDIASAK"        
      [13,] "sp|P07262|DHE4_YEAST" "GCIISETGITSEQVADISSAK"        
      [14,] "sp|P39708|DHE5_YEAST" "GGLCVDLK"                     
      [15,] "sp|P07262|DHE4_YEAST" "GGLCVDLK"                     
      [16,] "sp|P07262|DHE4_YEAST" "HIGQDTDVPAGDIGVGGR"           
      [17,] "sp|P39708|DHE5_YEAST" "ICYAFMR"                      
      [18,] "sp|P07262|DHE4_YEAST" "ICYAFMR"                      
      [19,] "sp|P07262|DHE4_YEAST" "IMINCFNECIDYAK"               
      [20,] "sp|P07262|DHE4_YEAST" "ITWTSER"                      
      [21,] "sp|P39708|DHE5_YEAST" "NSWEGVLTGK"                   
      [22,] "sp|P07262|DHE4_YEAST" "NSWEGVLTGK"                   
      [23,] "sp|P07262|DHE4_YEAST" "SEPEFQQAYEEVVSSLEDSTLFEQHPEYR"
      [24,] "sp|P39708|DHE5_YEAST" "SLEEIVDEYSTFSESK"             
      [25,] "sp|P07262|DHE4_YEAST" "SLEQIVNEYSTFSENK"             
      [26,] "sp|P07262|DHE4_YEAST" "STATGPSEAVWYGPPK"             
      [27,] "sp|P07262|DHE4_YEAST" "VDIALPCATQNEVSGEEAK"          
      [28,] "sp|P07262|DHE4_YEAST" "VIELGGTVVSLSDSK"              
      [29,] "sp|P39708|DHE5_YEAST" "VLPIVSVPER"                   
      [30,] "sp|P07262|DHE4_YEAST" "VLPIVSVPER"                   
      [31,] "sp|P39708|DHE5_YEAST" "VTISGSGNVAQYAALK"             
      [32,] "sp|P07262|DHE4_YEAST" "VTISGSGNVAQYAALK"             
      [33,] "sp|P07262|DHE4_YEAST" "VTWENDKGEQEVAQGYR"            
      [34,] "sp|P39708|DHE5_YEAST" "VTWENDNGEQEVAQGYR"            

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA  0.001710729  0.418949419 -0.281233039
       [6] -0.512467896 -0.003605408 -0.225957852  0.121464806 -0.239213508
      [11] -0.041718450 -0.247831253  0.015750668 -0.081089518  0.195036467
      [16] -0.002439029 -0.082496512 -0.023823956  0.189438722  0.208994714
      [21]  0.812429458  0.222805137  0.570468298 -0.288777121  0.206792322
      [26] -0.182535661 -0.335666498  0.377177558

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                          
       [1,] "sp|P09938|RIR2_YEAST" "AEWALR"                      
       [2,] "sp|P09938|RIR2_YEAST" "DIHDWNNR"                    
       [3,] "sp|P09938|RIR2_YEAST" "ELETLREENR"                  
       [4,] "sp|P09938|RIR2_YEAST" "ESEFLFNAIHTIPEIGEK"          
       [5,] "sp|P09938|RIR2_YEAST" "IVTEAVEIEQR"                 
       [6,] "sp|P09938|RIR2_YEAST" "LLVAFGNK"                    
       [7,] "sp|P09938|RIR2_YEAST" "NKPDPAIVEK"                  
       [8,] "sp|P09938|RIR2_YEAST" "TNFFEK"                      
       [9,] "sp|P09938|RIR2_YEAST" "TVLFPIK"                     
      [10,] "sp|P09938|RIR2_YEAST" "VENPFDFMENISLAGK"            
      [11,] "sp|P09938|RIR2_YEAST" "WIQDADALFGER"                
      [12,] "sp|P09938|RIR2_YEAST" "YFLDALPVALLGMNADLMNQYVEFVADR"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA  0.29906710  0.40978862 -0.59302233 -0.02610465 -0.37618753
       [7]  0.01345664 -0.05842363  0.22224993 -0.01225788  0.34024455  0.08348425
      [13] -0.21944505

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                          
       [1,] "sp|P00330|ADH1_YEAST" "ANELLINVK"                   
       [2,] "sp|P00331|ADH2_YEAST" "ANGTVVLVGLPAGAK"             
       [3,] "sp|P00330|ADH1_YEAST" "ATDGGAHGVINVSVSEAAIEASTR"    
       [4,] "sp|P00330|ADH1_YEAST" "CCSDVFNQVVK"                 
       [5,] "sp|P00331|ADH2_YEAST" "CSSDVFNHVVK"                 
       [6,] "sp|P07246|ADH3_YEAST" "DIPVPEPKPNEILINVK"           
       [7,] "sp|P00330|ADH1_YEAST" "DIVGAVLK"                    
       [8,] "sp|P00330|ADH1_YEAST" "EALDFFAR"                    
       [9,] "sp|P00331|ADH2_YEAST" "EALDFFAR"                    
      [10,] "sp|P38113|ADH5_YEAST" "EALDFFAR"                    
      [11,] "sp|P07246|ADH3_YEAST" "EALDFFSR"                    
      [12,] "sp|P00330|ADH1_YEAST" "EKDIVGAVLK"                  
      [13,] "sp|P07246|ADH3_YEAST" "GVIFYENK"                    
      [14,] "sp|P00330|ADH1_YEAST" "GVIFYESHGK"                  
      [15,] "sp|P00330|ADH1_YEAST" "IGDYAGIK"                    
      [16,] "sp|P00331|ADH2_YEAST" "IGDYAGIK"                    
      [17,] "sp|P07246|ADH3_YEAST" "IQQGTDLAEVAPILCAGVTVYK"      
      [18,] "sp|P07246|ADH3_YEAST" "IVGLSELPK"                   
      [19,] "sp|P00330|ADH1_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [20,] "sp|P00331|ADH2_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [21,] "sp|P07246|ADH3_YEAST" "NMVSDIQEATK"                 
      [22,] "sp|P00330|ADH1_YEAST" "SANLMAGHWVAISGAAGGLGSLAVQYAK"
      [23,] "sp|P00330|ADH1_YEAST" "SIGGEVFIDFTK"                
      [24,] "sp|P00330|ADH1_YEAST" "SIPETQK"                     
      [25,] "sp|P00331|ADH2_YEAST" "SIPETQK"                     
      [26,] "sp|P00330|ADH1_YEAST" "SISIVGSYVGNR"                
      [27,] "sp|P00331|ADH2_YEAST" "SISIVGSYVGNR"                
      [28,] "sp|P00330|ADH1_YEAST" "VLGIDGGEGK"                  
      [29,] "sp|P00330|ADH1_YEAST" "VLGIDGGEGKEELFR"             
      [30,] "sp|P00331|ADH2_YEAST" "VVGLSSLPEIYEK"               
      [31,] "sp|P00330|ADH1_YEAST" "VVGLSTLPEIYEK"               
      [32,] "sp|P00330|ADH1_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       
      [33,] "sp|P07246|ADH3_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA           NA           NA  0.047244275
       [6]  0.156264076 -0.308834358 -0.122943668  0.052899789 -0.263012183
      [11]  0.063686932  0.027271420  0.381499448 -0.058657242  0.089598227
      [16] -0.902526234  0.193407384  0.315243882  0.274948419  0.007575608
      [21]  0.016009791 -0.607872347  0.059782914 -0.533626424  0.045233495
      [26]  0.142663090 -0.318964500  0.245981868  0.095958510  0.220024851

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 10

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA  0.21847490  2.38974845 -0.92220549  0.06547634
       [7]  0.01649847  0.64962728  1.57342573  0.18904660 -0.37392833 -0.42862844

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 26

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA  0.001710729  0.418949419 -0.281233039
       [6] -0.512467896 -0.003605408 -0.225957852  0.121464806 -0.239213508
      [11] -0.041718450 -0.247831253  0.015750668 -0.081089518  0.195036467
      [16] -0.002439029 -0.082496512 -0.023823956  0.189438722  0.208994714
      [21]  0.812429458  0.222805137  0.570468298 -0.288777121  0.206792322
      [26] -0.182535661 -0.335666498  0.377177558

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 12

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA  0.29906710  0.40978862 -0.59302233 -0.02610465 -0.37618753
       [7]  0.01345664 -0.05842363  0.22224993 -0.01225788  0.34024455  0.08348425
      [13] -0.21944505

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 26

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA           NA           NA  0.047244275
       [6]  0.156264076 -0.308834358 -0.122943668  0.052899789 -0.263012183
      [11]  0.063686932  0.027271420  0.381499448 -0.058657242  0.089598227
      [16] -0.902526234  0.193407384  0.315243882  0.274948419  0.007575608
      [21]  0.016009791 -0.607872347  0.059782914 -0.533626424  0.045233495
      [26]  0.142663090 -0.318964500  0.245981868  0.095958510  0.220024851

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

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
       [1]           NA           NA  0.287692950  2.264404689 -0.286033804
       [6] -0.079327006  0.009683422  1.006707013  1.647448282  0.232567540
      [11] -0.335984406 -0.435527101

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                           
       [1,] "sp|P39708|DHE5_YEAST" "AANLGGVAVSGLEMAQNSQK"         
       [2,] "sp|P07262|DHE4_YEAST" "AANLGGVAVSGLEMAQNSQR"         
       [3,] "sp|P39708|DHE5_YEAST" "FHPSVNLSILK"                  
       [4,] "sp|P07262|DHE4_YEAST" "FHPSVNLSILK"                  
       [5,] "sp|P07262|DHE4_YEAST" "FIAEGSNMGSTPEAIAVFETAR"       
       [6,] "sp|P39708|DHE5_YEAST" "FLGFEQIFK"                    
       [7,] "sp|P07262|DHE4_YEAST" "FLGFEQIFK"                    
       [8,] "sp|P07262|DHE4_YEAST" "GANIASFIK"                    
       [9,] "sp|P39708|DHE5_YEAST" "GANIASFVMVADAMLDQGDVF"        
      [10,] "sp|P39708|DHE5_YEAST" "GCIISETGITSEQIHDIASAK"        
      [11,] "sp|P07262|DHE4_YEAST" "GCIISETGITSEQVADISSAK"        
      [12,] "sp|P39708|DHE5_YEAST" "GGLCVDLK"                     
      [13,] "sp|P07262|DHE4_YEAST" "GGLCVDLK"                     
      [14,] "sp|P07262|DHE4_YEAST" "HIGQDTDVPAGDIGVGGR"           
      [15,] "sp|P39708|DHE5_YEAST" "ICYAFMR"                      
      [16,] "sp|P07262|DHE4_YEAST" "ICYAFMR"                      
      [17,] "sp|P07262|DHE4_YEAST" "IMINCFNECIDYAK"               
      [18,] "sp|P07262|DHE4_YEAST" "ITWTSER"                      
      [19,] "sp|P39708|DHE5_YEAST" "NSWEGVLTGK"                   
      [20,] "sp|P07262|DHE4_YEAST" "NSWEGVLTGK"                   
      [21,] "sp|P07262|DHE4_YEAST" "SEPEFQQAYEEVVSSLEDSTLFEQHPEYR"
      [22,] "sp|P39708|DHE5_YEAST" "SLEEIVDEYSTFSESK"             
      [23,] "sp|P07262|DHE4_YEAST" "SLEQIVNEYSTFSENK"             
      [24,] "sp|P07262|DHE4_YEAST" "STATGPSEAVWYGPPK"             
      [25,] "sp|P07262|DHE4_YEAST" "VDIALPCATQNEVSGEEAK"          
      [26,] "sp|P07262|DHE4_YEAST" "VIELGGTVVSLSDSK"              
      [27,] "sp|P39708|DHE5_YEAST" "VLPIVSVPER"                   
      [28,] "sp|P07262|DHE4_YEAST" "VLPIVSVPER"                   
      [29,] "sp|P07262|DHE4_YEAST" "VTWENDKGEQEVAQGYR"            
      [30,] "sp|P39708|DHE5_YEAST" "VTWENDNGEQEVAQGYR"            

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA  0.133767606 -0.102232106 -0.707239921
       [6]  0.067674521 -0.313992414  0.264007794 -0.485394528 -0.108570991
      [11] -0.252855142 -0.162971396 -0.020649276  0.237542543 -0.090255292
      [16]  0.008389828 -0.078087345  0.073598793  0.254094730  0.640924076
      [21]  0.225849234  0.724753721 -0.229717380  0.388979230 -0.378671752
      [26]  0.391706817

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                          
       [1,] "sp|P09938|RIR2_YEAST" "AEWALR"                      
       [2,] "sp|P09938|RIR2_YEAST" "ELETLREENR"                  
       [3,] "sp|P09938|RIR2_YEAST" "ESEFLFNAIHTIPEIGEK"          
       [4,] "sp|P09938|RIR2_YEAST" "IVTEAVEIEQR"                 
       [5,] "sp|P09938|RIR2_YEAST" "LLVAFGNK"                    
       [6,] "sp|P09938|RIR2_YEAST" "NKPDPAIVEK"                  
       [7,] "sp|P09938|RIR2_YEAST" "TNFFEK"                      
       [8,] "sp|P09938|RIR2_YEAST" "TVLFPIK"                     
       [9,] "sp|P09938|RIR2_YEAST" "VENPFDFMENISLAGK"            
      [10,] "sp|P09938|RIR2_YEAST" "WIQDADALFGER"                
      [11,] "sp|P09938|RIR2_YEAST" "YFLDALPVALLGMNADLMNQYVEFVADR"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA  0.28012701 -0.40490051  0.08369477 -0.46353613  0.25279671
       [7] -0.01930214  0.26473825 -0.03570867  0.16171324  0.03883347 -0.18292606

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                          
       [1,] "sp|P00330|ADH1_YEAST" "ANELLINVK"                   
       [2,] "sp|P00331|ADH2_YEAST" "ANGTVVLVGLPAGAK"             
       [3,] "sp|P00330|ADH1_YEAST" "ATDGGAHGVINVSVSEAAIEASTR"    
       [4,] "sp|P00330|ADH1_YEAST" "CCSDVFNQVVK"                 
       [5,] "sp|P00331|ADH2_YEAST" "CSSDVFNHVVK"                 
       [6,] "sp|P07246|ADH3_YEAST" "DIPVPEPKPNEILINVK"           
       [7,] "sp|P00330|ADH1_YEAST" "DIVGAVLK"                    
       [8,] "sp|P00330|ADH1_YEAST" "EALDFFAR"                    
       [9,] "sp|P00331|ADH2_YEAST" "EALDFFAR"                    
      [10,] "sp|P38113|ADH5_YEAST" "EALDFFAR"                    
      [11,] "sp|P07246|ADH3_YEAST" "EALDFFSR"                    
      [12,] "sp|P00330|ADH1_YEAST" "EKDIVGAVLK"                  
      [13,] "sp|P07246|ADH3_YEAST" "GVIFYENK"                    
      [14,] "sp|P00330|ADH1_YEAST" "GVIFYESHGK"                  
      [15,] "sp|P00330|ADH1_YEAST" "IGDYAGIK"                    
      [16,] "sp|P00331|ADH2_YEAST" "IGDYAGIK"                    
      [17,] "sp|P07246|ADH3_YEAST" "IQQGTDLAEVAPILCAGVTVYK"      
      [18,] "sp|P07246|ADH3_YEAST" "IVGLSELPK"                   
      [19,] "sp|P00330|ADH1_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [20,] "sp|P00331|ADH2_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [21,] "sp|P07246|ADH3_YEAST" "NMVSDIQEATK"                 
      [22,] "sp|P00330|ADH1_YEAST" "SANLMAGHWVAISGAAGGLGSLAVQYAK"
      [23,] "sp|P00330|ADH1_YEAST" "SIGGEVFIDFTK"                
      [24,] "sp|P00330|ADH1_YEAST" "SIPETQK"                     
      [25,] "sp|P00331|ADH2_YEAST" "SIPETQK"                     
      [26,] "sp|P00330|ADH1_YEAST" "SISIVGSYVGNR"                
      [27,] "sp|P00331|ADH2_YEAST" "SISIVGSYVGNR"                
      [28,] "sp|P00330|ADH1_YEAST" "VLGIDGGEGK"                  
      [29,] "sp|P00330|ADH1_YEAST" "VLGIDGGEGKEELFR"             
      [30,] "sp|P00331|ADH2_YEAST" "VVGLSSLPEIYEK"               
      [31,] "sp|P00330|ADH1_YEAST" "VVGLSTLPEIYEK"               
      [32,] "sp|P00330|ADH1_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       
      [33,] "sp|P07246|ADH3_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA          NA          NA -0.04603679  0.24473348
       [7] -0.06265124 -0.27922238 -0.54851506  0.25606028 -0.06370699  0.04999654
      [13]  0.46765717 -0.12193354 -0.02614743 -0.96783305  0.06120803  0.26019685
      [19]  0.36389932  0.25941473 -0.01404349 -0.38089762  0.17378518 -0.85100911
      [25]  0.08933018 -0.07592864 -0.33853470  0.16662052  0.09253790  0.18161783

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 10

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA  0.287692950  2.264404689 -0.286033804
       [6] -0.079327006  0.009683422  1.006707013  1.647448282  0.232567540
      [11] -0.335984406 -0.435527101

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 24

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA  0.133767606 -0.102232106 -0.707239921
       [6]  0.067674521 -0.313992414  0.264007794 -0.485394528 -0.108570991
      [11] -0.252855142 -0.162971396 -0.020649276  0.237542543 -0.090255292
      [16]  0.008389828 -0.078087345  0.073598793  0.254094730  0.640924076
      [21]  0.225849234  0.724753721 -0.229717380  0.388979230 -0.378671752
      [26]  0.391706817

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 11

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA  0.28012701 -0.40490051  0.08369477 -0.46353613  0.25279671
       [7] -0.01930214  0.26473825 -0.03570867  0.16171324  0.03883347 -0.18292606

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 26

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA          NA          NA -0.04603679  0.24473348
       [7] -0.06265124 -0.27922238 -0.54851506  0.25606028 -0.06370699  0.04999654
      [13]  0.46765717 -0.12193354 -0.02614743 -0.96783305  0.06120803  0.26019685
      [19]  0.36389932  0.25941473 -0.01404349 -0.38089762  0.17378518 -0.85100911
      [25]  0.08933018 -0.07592864 -0.33853470  0.16662052  0.09253790  0.18161783

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

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
       [1]          NA          NA  0.12101599  3.01458538 -0.32998840 -0.14905643
       [7] -0.03189065 -0.03171039  1.90939596  0.24880656 -0.29682901 -0.54509853

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                           
       [1,] "sp|P39708|DHE5_YEAST" "AANLGGVAVSGLEMAQNSQK"         
       [2,] "sp|P07262|DHE4_YEAST" "AANLGGVAVSGLEMAQNSQR"         
       [3,] "sp|P39708|DHE5_YEAST" "FHPSVNLSILK"                  
       [4,] "sp|P07262|DHE4_YEAST" "FHPSVNLSILK"                  
       [5,] "sp|P07262|DHE4_YEAST" "FIAEGSNMGSTPEAIAVFETAR"       
       [6,] "sp|P39708|DHE5_YEAST" "FLGFEQIFK"                    
       [7,] "sp|P07262|DHE4_YEAST" "FLGFEQIFK"                    
       [8,] "sp|P07262|DHE4_YEAST" "GANIASFIK"                    
       [9,] "sp|P39708|DHE5_YEAST" "GANIASFVMVADAMLDQGDVF"        
      [10,] "sp|P39708|DHE5_YEAST" "GCIISETGITSEQIHDIASAK"        
      [11,] "sp|P07262|DHE4_YEAST" "GCIISETGITSEQVADISSAK"        
      [12,] "sp|P39708|DHE5_YEAST" "GGLCVDLK"                     
      [13,] "sp|P07262|DHE4_YEAST" "GGLCVDLK"                     
      [14,] "sp|P07262|DHE4_YEAST" "HIGQDTDVPAGDIGVGGR"           
      [15,] "sp|P39708|DHE5_YEAST" "ICYAFMR"                      
      [16,] "sp|P07262|DHE4_YEAST" "ICYAFMR"                      
      [17,] "sp|P07262|DHE4_YEAST" "IMINCFNECIDYAK"               
      [18,] "sp|P07262|DHE4_YEAST" "ITWTSER"                      
      [19,] "sp|P39708|DHE5_YEAST" "NSWEGVLTGK"                   
      [20,] "sp|P07262|DHE4_YEAST" "NSWEGVLTGK"                   
      [21,] "sp|P07262|DHE4_YEAST" "SEPEFQQAYEEVVSSLEDSTLFEQHPEYR"
      [22,] "sp|P39708|DHE5_YEAST" "SLEEIVDEYSTFSESK"             
      [23,] "sp|P07262|DHE4_YEAST" "SLEQIVNEYSTFSENK"             
      [24,] "sp|P07262|DHE4_YEAST" "STATGPSEAVWYGPPK"             
      [25,] "sp|P07262|DHE4_YEAST" "VDIALPCATQNEVSGEEAK"          
      [26,] "sp|P07262|DHE4_YEAST" "VIELGGTVVSLSDSK"              
      [27,] "sp|P39708|DHE5_YEAST" "VLPIVSVPER"                   
      [28,] "sp|P07262|DHE4_YEAST" "VLPIVSVPER"                   
      [29,] "sp|P39708|DHE5_YEAST" "VTISGSGNVAQYAALK"             
      [30,] "sp|P07262|DHE4_YEAST" "VTISGSGNVAQYAALK"             
      [31,] "sp|P07262|DHE4_YEAST" "VTWENDKGEQEVAQGYR"            
      [32,] "sp|P39708|DHE5_YEAST" "VTWENDNGEQEVAQGYR"            

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA  0.108869474  0.266656696 -0.423961204
       [6] -0.009416542 -0.366922769  0.077044870 -0.762952515 -0.179557111
      [11] -0.326441969 -0.161826295 -0.152779282 -0.144939830 -0.059443388
      [16] -0.128216106 -0.149743098  0.033196629  0.317197786  0.797830589
      [21]  0.394707012  0.798598095 -0.264435785  0.453109316  0.615788613
      [26] -0.539799757  0.522376098

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                          
       [1,] "sp|P09938|RIR2_YEAST" "AEWALR"                      
       [2,] "sp|P09938|RIR2_YEAST" "DEGLHTDFACLLFAHLK"           
       [3,] "sp|P09938|RIR2_YEAST" "ESEFLFNAIHTIPEIGEK"          
       [4,] "sp|P09938|RIR2_YEAST" "IVTEAVEIEQR"                 
       [5,] "sp|P09938|RIR2_YEAST" "LLVAFGNK"                    
       [6,] "sp|P09938|RIR2_YEAST" "NKPDPAIVEK"                  
       [7,] "sp|P09938|RIR2_YEAST" "TNFFEK"                      
       [8,] "sp|P09938|RIR2_YEAST" "TVLFPIK"                     
       [9,] "sp|P09938|RIR2_YEAST" "VENPFDFMENISLAGK"            
      [10,] "sp|P09938|RIR2_YEAST" "WIQDADALFGER"                
      [11,] "sp|P09938|RIR2_YEAST" "YFLDALPVALLGMNADLMNQYVEFVADR"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA  0.19122515 -0.52578414  0.10908731 -0.39811599  0.55840120
       [7] -0.29849697  0.27061181  0.09033515  0.17864350  0.42955369 -0.18767689

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
       [1]          NA          NA          NA          NA -0.05716547  0.08877794
       [7] -0.10061788 -0.17422051  0.44044457 -0.09454628  0.06518158  0.41923160
      [13] -0.18761462 -0.06393003 -1.03477697  0.24399061  0.09832052  0.31545243
      [19]  0.19333349  0.11567471 -0.58879066  0.20612611 -0.39382336  0.09875934
      [25]  1.04391822 -0.42545315  0.33826907  0.09983684  0.68865357

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 10

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA  0.12101599  3.01458538 -0.32998840 -0.14905643
       [7] -0.03189065 -0.03171039  1.90939596  0.24880656 -0.29682901 -0.54509853

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
       [1]           NA           NA  0.108869474  0.266656696 -0.423961204
       [6] -0.009416542 -0.366922769  0.077044870 -0.762952515 -0.179557111
      [11] -0.326441969 -0.161826295 -0.152779282 -0.144939830 -0.059443388
      [16] -0.128216106 -0.149743098  0.033196629  0.317197786  0.797830589
      [21]  0.394707012  0.798598095 -0.264435785  0.453109316  0.615788613
      [26] -0.539799757  0.522376098

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

---

    Code
      sum(!igraph::V(graphsImp[[i]][[j]])$type)
    Output
      [1] 11

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA  0.19122515 -0.52578414  0.10908731 -0.39811599  0.55840120
       [7] -0.29849697  0.27061181  0.09033515  0.17864350  0.42955369 -0.18767689

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
       [1]          NA          NA          NA          NA -0.05716547  0.08877794
       [7] -0.10061788 -0.17422051  0.44044457 -0.09454628  0.06518158  0.41923160
      [13] -0.18761462 -0.06393003 -1.03477697  0.24399061  0.09832052  0.31545243
      [19]  0.19333349  0.11567471 -0.58879066  0.20612611 -0.39382336  0.09875934
      [25]  1.04391822 -0.42545315  0.33826907  0.09983684  0.68865357

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      NULL

