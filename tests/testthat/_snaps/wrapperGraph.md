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
       [1]          NA          NA  0.20055863  2.60772144 -0.10455462  0.57623097
       [7] -0.02534580  1.00071495  0.91778084  0.08702881 -0.63655573

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
       [1]          NA          NA  0.12523316 -0.14016680 -0.11745684 -0.72738638
       [7]  0.10512926 -0.22006861  0.12462445 -0.22891262 -0.03882496 -0.66511747
      [13] -0.11369859  0.39197902 -0.03849029 -0.61229730  0.24022045 -0.30390075
      [19]  0.37879733  0.53215973  0.15132580  0.58982476 -0.24773651  0.12652204
      [25]  0.72859391 -0.26178118  0.16356172

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
       [1]          NA  0.33734292 -0.80011053  0.61113343 -0.18172860  0.15185502
       [7] -0.42718020  0.11935644 -0.03393225 -0.10305954  0.02370556  0.14843254
      [13]  0.27604730  0.17149957 -0.34449242

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
       [1]           NA           NA           NA           NA  0.034291942
       [6]  0.202500716 -0.315618274 -0.118970527  0.162209921  0.282490167
      [11]  0.002504345  0.024930690  0.283225570 -0.228069875  0.048889130
      [16] -0.719635602  0.165156750  0.202055936  0.302234371  0.032762446
      [21] -0.601152274 -0.714167205  0.316199116  0.052617110  0.073263795
      [26]  0.217881661 -0.207179833  0.230756591  0.077885092  0.120472806

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
       [1]           NA           NA  0.220744867  2.354046114 -0.805105278
       [6]  0.034852845 -0.009964474  0.666194811  1.560851490  0.199090597
      [11] -0.404801795 -0.311993423

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
       [1]          NA          NA  0.04155965  0.38966718 -0.30002825 -0.53104674
       [7] -0.02547976 -0.24698979  0.11132551 -0.23369879 -0.02063805 -0.28550356
      [13]  0.12345178 -0.10655481  0.27123025  0.10681787 -0.03129128  0.08053199
      [19]  0.19765335  0.23783725  0.79609537  0.20284198  0.55436279 -0.19208566
      [25]  0.18530586 -0.31650176 -0.35923144  0.32109131

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
       [1]          NA  0.27738920  0.36285289 -0.61825026 -0.05141660 -0.40507114
       [7]  0.07984209 -0.06551464  0.20715715  0.10296698  0.37670346  0.10104503
      [13] -0.24943864

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
       [1]           NA           NA           NA           NA  0.056157559
       [6]  0.131311411 -0.300902827 -0.134982549  0.017984400 -0.233384153
      [11]  0.061608147  0.038493401  0.371295077 -0.083471067  0.061838751
      [16] -0.898266537  0.195154050  0.324611138  0.261309738  0.012911132
      [21] -0.002535603 -0.636321043  0.066915559 -0.564194609  0.061057662
      [26]  0.205973487 -0.302934899  0.224149400  0.113110891  0.223350052

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
       [1]          NA          NA  0.29288178  2.23673999 -0.15770229 -0.11937266
       [7] -0.03300617  1.02949695  1.62133083  0.24408006 -0.37108360 -0.31119211

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
       [1]          NA          NA  0.17975584 -0.14223025 -0.72927442  0.03420529
       [7] -0.33485315  0.26191968 -0.47708036 -0.07927414 -0.32666845 -0.04982611
      [13] -0.06084031  0.32681754  0.02802284  0.06774371  0.03521061  0.08335819
      [19]  0.29062465  0.63090940  0.20627259  0.71244392 -0.12576884  0.36035698
      [25] -0.41092183  0.28362848

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
       [1]          NA  0.24475589 -0.45395714  0.04046221 -0.48536629  0.33265624
       [7] -0.01902505  0.24375845  0.08966783  0.20662622  0.06419420 -0.22022447

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
       [1]           NA           NA           NA           NA -0.033168766
       [6]  0.201969222 -0.054314513 -0.284761467 -0.607336035  0.284044435
      [11] -0.069635987  0.065758004  0.463134245 -0.160286467 -0.056277472
      [16] -0.961210805  0.063567391  0.276858164  0.357511337  0.267320402
      [21] -0.036806731 -0.423354873  0.181973587 -0.886279392  0.110121805
      [26] -0.003981729 -0.317285140  0.130767589  0.115346289  0.188783820

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
       [1]          NA          NA  0.12372800  2.99048809 -0.19995070 -0.17330655
       [7] -0.06909708 -0.01290987  1.87704339  0.26561986 -0.31996035 -0.43455735

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
       [1]          NA          NA  0.14570322  0.23510602 -0.44903269 -0.05327657
       [7] -0.40527897  0.07508703 -0.75389543 -0.15949689 -0.39140950 -0.04579429
      [13] -0.19478058 -0.06082763  0.06011639 -0.06699893 -0.03368252  0.04512118
      [19]  0.34648305  0.77390590  0.35745521  0.77019016 -0.15651843  0.41137986
      [25]  0.69027188 -0.58415670  0.43068296

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
       [1]         NA  0.1552400 -0.5838083  0.0695772 -0.4174585  0.6295286
       [7] -0.2988141  0.2464730  0.2187936  0.2222274  0.4547324 -0.2121912

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
       [1]          NA          NA          NA          NA -0.04076186  0.04793029
       [7] -0.08517390 -0.19163029  0.46810543 -0.11109183  0.08491465  0.40895239
      [13] -0.23339927 -0.09617435 -1.02396070  0.24925780  0.10723328  0.30957137
      [19]  0.20569944  0.08773296 -0.61427416  0.21718992 -0.40769424  0.12418496
      [25]  1.11515574 -0.40064889  0.29540792  0.12463580  0.69749011

