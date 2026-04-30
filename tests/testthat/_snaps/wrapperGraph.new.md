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
            [,1]                  
       [1,] "sp|P07262|DHE4_YEAST"
       [2,] "sp|P07262|DHE4_YEAST"
       [3,] "sp|P07262|DHE4_YEAST"
       [4,] "sp|P07262|DHE4_YEAST"
       [5,] "sp|P07262|DHE4_YEAST"
       [6,] "sp|P07262|DHE4_YEAST"
       [7,] "sp|P07262|DHE4_YEAST"
       [8,] "sp|P07262|DHE4_YEAST"
       [9,] "sp|P07262|DHE4_YEAST"
      [10,] "sp|P07262|DHE4_YEAST"
      [11,] "sp|P07262|DHE4_YEAST"
      [12,] "sp|P07262|DHE4_YEAST"
      [13,] "sp|P07262|DHE4_YEAST"
      [14,] "sp|P07262|DHE4_YEAST"
      [15,] "sp|P07262|DHE4_YEAST"
      [16,] "sp|P07262|DHE4_YEAST"
      [17,] "sp|P07262|DHE4_YEAST"
      [18,] "sp|P07262|DHE4_YEAST"
      [19,] "sp|P39708|DHE5_YEAST"
      [20,] "sp|P39708|DHE5_YEAST"
      [21,] "sp|P39708|DHE5_YEAST"
      [22,] "sp|P39708|DHE5_YEAST"
      [23,] "sp|P39708|DHE5_YEAST"
      [24,] "sp|P39708|DHE5_YEAST"
      [25,] "sp|P39708|DHE5_YEAST"
      [26,] "sp|P39708|DHE5_YEAST"
      [27,] "sp|P39708|DHE5_YEAST"
      [28,] "sp|P39708|DHE5_YEAST"
      [29,] "sp|P39708|DHE5_YEAST"
      [30,] "sp|P39708|DHE5_YEAST"
      [31,] "sp|P39708|DHE5_YEAST"
      [32,] "sp|P39708|DHE5_YEAST"
      [33,] "sp|P39708|DHE5_YEAST"
      [34,] "sp|P39708|DHE5_YEAST"
      [35,] "sp|P39708|DHE5_YEAST"
            [,2]                                                
       [1,] "AANLGGVAVSGLEMAQNSQR"                              
       [2,] "ALVAQGVK"                                          
       [3,] "ALVAQGVKFIAEGSNMGSTPEAIAVFETARSTATGPSEAVWYGPPK"    
       [4,] "GGLRFHPSVNLSILK"                                   
       [5,] "GPYKGGLR"                                          
       [6,] "GPYKGGLRFHPSVNLSILK"                               
       [7,] "KVLPIVSVPERIIQFR"                                  
       [8,] "NSWEGVLTGKGLNWGGSLIRPEATGYGLVYYTQAMIDYATNGK"       
       [9,] "RVTISGSGNVAQYAALK"                                 
      [10,] "SLEQIVNEYSTFSENKVQYIAGARPWTHVQKVDIALPCATQNEVSGEEAK"
      [11,] "SNNEIRRICYAFMR"                                    
      [12,] "VDQELKR"                                           
      [13,] "VQYIAGARPWTHVQK"                                   
      [14,] "VQYIAGARPWTHVQKVDIALPCATQNEVSGEEAKALVAQGVK"        
      [15,] "VQYNSAK"                                           
      [16,] "VQYNSAKGPYK"                                       
      [17,] "VTISGSGNVAQYAALK"                                  
      [18,] "VTWENDKGEQEVAQGYR"                                 
      [19,] "DTDVPAGDIGVGGR"                                    
      [20,] "FKSLEEIVDEYSTFSESKMK"                              
      [21,] "GGLRFHPSVNLSILK"                                   
      [22,] "GPYKGGLR"                                          
      [23,] "GPYKGGLRFHPSVNLSILK"                               
      [24,] "HIGKDTDVPAGDIGVGGREIGYLFGAYR"                      
      [25,] "IFEKFPQYK"                                         
      [26,] "KIMINCFNDCIQAAQEYSTEK"                             
      [27,] "KVLPIVSVPERIIQFR"                                  
      [28,] "MKYVAGARPWTHVSNVDIALPCATQNEVSGDEAKALVASGVK"        
      [29,] "RVTISGSGNVAQYAALK"                                 
      [30,] "SDNEIRR"                                           
      [31,] "SLEEIVDEYSTFSESK"                                  
      [32,] "SYKNSWEGVLTGKGLNWGGSLIRPEATGFGLVYYTQAMIDYATNGK"    
      [33,] "TSEPEFQQAYDEIVSSVEDSK"                             
      [34,] "VTISGSGNVAQYAALK"                                  
      [35,] "VTWENDNGEQEVAQGYRVQFNSAK"                          

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA  0.024154645  0.045874025 -0.082146399
       [6]  0.001562107  0.020865276  0.109482415 -0.056102585 -0.001518065
      [11] -0.002770665  0.032973967 -0.018048014  0.031696259  0.051000556
      [16]  0.102364789  0.048377111  0.037552670 -0.008033538 -0.021723577
      [21]  0.105190747  0.082825205 -0.079188015 -0.077989000 -0.009707192
      [26]  0.063765156 -0.021042240  0.163287032 -0.012062660 -0.102809834
      [31] -0.031454281

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                                          [,2]                 
       [1,] "sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST" "AEAEAEKK"           
       [2,] "sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST" "AISKNLPILKNHFR"     
       [3,] "sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST" "EKAEAEAEK"          
       [4,] "sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST" "GIREKR"             
       [5,] "sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST" "HWQERVKVHFDQAGK"    
       [6,] "sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST" "IAPRPLDLLRPVVRAPTVK"
       [7,] "sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST" "NLPILKNHFR"         
       [8,] "sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST" "NQEIFDANVQRLK"      
       [9,] "sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST" "TIGIAVDHR"          
      [10,] "sp|P40212|RL13B_YEAST;sp|Q12690|RL13A_YEAST" "TLRLARSEK"          

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA -0.041732672  0.004575050 -0.165466168 -0.144230445
       [6] -0.038326070  0.115425374 -0.068608730 -0.008191886  0.092838731
      [11] -0.038088544

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                                     
       [1,] "sp|P09938|RIR2_YEAST" "DEGLHTDFACLLFAHLK"                      
       [2,] "sp|P09938|RIR2_YEAST" "DIHDWNNRMNENER"                         
       [3,] "sp|P09938|RIR2_YEAST" "DIHDWNNRMNENERFFISR"                    
       [4,] "sp|P09938|RIR2_YEAST" "ERTVLFPIKYHEIWQAYK"                     
       [5,] "sp|P09938|RIR2_YEAST" "IVTEAVEIEQRYFLDALPVALLGMNADLMNQYVEFVADR"
       [6,] "sp|P09938|RIR2_YEAST" "LLVAFGNK"                               
       [7,] "sp|P09938|RIR2_YEAST" "LLVAFGNKK"                              
       [8,] "sp|P09938|RIR2_YEAST" "LSKDAENHKAYLK"                          
       [9,] "sp|P09938|RIR2_YEAST" "SDMLKEKLSK"                             
      [10,] "sp|P09938|RIR2_YEAST" "WIQDADALFGERLVAFASIEGVFFSGSFASIFWLK"    

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA  0.06716363  0.03344063  0.11232045 -0.08385854 -0.07405091
       [7]  0.07122601 -0.00350411  0.09586286 -0.07074105 -0.01626986

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                    [,2]                         
       [1,] "sp|P40212|RL13B_YEAST" "AAKIAPRPLDLLRPVVR"          
       [2,] "sp|Q12690|RL13A_YEAST" "AAKIAPRPLDLLRPVVR"          
       [3,] "sp|P40212|RL13B_YEAST" "AEAEAEK"                    
       [4,] "sp|Q12690|RL13A_YEAST" "AEAEAEK"                    
       [5,] "sp|P40212|RL13B_YEAST" "AISKNLPILKNHFR"             
       [6,] "sp|Q12690|RL13A_YEAST" "AISKNLPILKNHFR"             
       [7,] "sp|P40212|RL13B_YEAST" "AVQDNGESAFR"                
       [8,] "sp|Q12690|RL13A_YEAST" "AVQDNGESAFR"                
       [9,] "sp|P40212|RL13B_YEAST" "EKAEAEAEK"                  
      [10,] "sp|Q12690|RL13A_YEAST" "EKAEAEAEK"                  
      [11,] "sp|P40212|RL13B_YEAST" "EYQSKIIVFPR"                
      [12,] "sp|Q12690|RL13A_YEAST" "EYQSKIIVFPR"                
      [13,] "sp|P40212|RL13B_YEAST" "GFTLAEVKAAGLTAAYARTIGIAVDHR"
      [14,] "sp|Q12690|RL13A_YEAST" "GFTLAEVKAAGLTAAYARTIGIAVDHR"
      [15,] "sp|P40212|RL13B_YEAST" "IAPRPLDLLRPVVRAPTVKYNR"     
      [16,] "sp|Q12690|RL13A_YEAST" "IAPRPLDLLRPVVRAPTVKYNR"     
      [17,] "sp|P40212|RL13B_YEAST" "MAISKNLPILKNHFR"            
      [18,] "sp|Q12690|RL13A_YEAST" "MAISKNLPILKNHFR"            
      [19,] "sp|P40212|RL13B_YEAST" "NLPILKNHFR"                 
      [20,] "sp|Q12690|RL13A_YEAST" "NLPILKNHFR"                 
      [21,] "sp|P40212|RL13B_YEAST" "RNARAAR"                    
      [22,] "sp|Q12690|RL13A_YEAST" "RNARATR"                    
      [23,] "sp|P40212|RL13B_YEAST" "TIGIAVDHR"                  
      [24,] "sp|Q12690|RL13A_YEAST" "TIGIAVDHR"                  
      [25,] "sp|P40212|RL13B_YEAST" "TIGIAVDHRR"                 
      [26,] "sp|Q12690|RL13A_YEAST" "TIGIAVDHRR"                 
      [27,] "sp|P40212|RL13B_YEAST" "TLRLARSEK"                  
      [28,] "sp|Q12690|RL13A_YEAST" "TLRLARSEK"                  

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA          NA -0.02244962 -0.01753529 -0.02966866 -0.02245139
       [7] -0.10020222  0.04367164  0.07422895 -0.02676368  0.07575796 -0.05223501
      [13] -0.09438190 -0.12859991 -0.01510093 -0.13619467 -0.01624250

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                  
       [1,] "sp|P09938|RIR2_YEAST"
       [2,] "sp|P09938|RIR2_YEAST"
       [3,] "sp|P09938|RIR2_YEAST"
       [4,] "sp|P09938|RIR2_YEAST"
       [5,] "sp|P09938|RIR2_YEAST"
       [6,] "sp|P09938|RIR2_YEAST"
       [7,] "sp|P09938|RIR2_YEAST"
       [8,] "sp|P09938|RIR2_YEAST"
       [9,] "sp|P09938|RIR2_YEAST"
      [10,] "sp|P09938|RIR2_YEAST"
      [11,] "sp|P09938|RIR2_YEAST"
      [12,] "sp|P09938|RIR2_YEAST"
      [13,] "sp|P09938|RIR2_YEAST"
      [14,] "sp|P09938|RIR2_YEAST"
      [15,] "sp|P09938|RIR2_YEAST"
      [16,] "sp|P09938|RIR2_YEAST"
            [,2]                                               
       [1,] "DIHDWNNRMNENERFFISR"                              
       [2,] "ELETLR"                                           
       [3,] "ERTVLFPIK"                                        
       [4,] "ERTVLFPIKYHEIWQAYK"                               
       [5,] "FFISRVLAFFAASDGIVNENLVENFSTEVQIPEAK"              
       [6,] "IVTEAVEIEQRYFLDALPVALLGMNADLMNQYVEFVADRLLVAFGNK"  
       [7,] "LLVAFGNKK"                                        
       [8,] "LSKDAENHKAYLK"                                    
       [9,] "NKPDPAIVEKIVTEAVEIEQRYFLDALPVALLGMNADLMNQYVEFVADR"
      [10,] "PKETPSKAAADALSDLEIK"                              
      [11,] "RGMMPGLTFSNELICR"                                 
      [12,] "SDMLKEKLSK"                                       
      [13,] "STKQEAGAFTFNEDF"                                  
      [14,] "TNFFEKRVSDYQK"                                    
      [15,] "VLAFFAASDGIVNENLVENFSTEVQIPEAK"                   
      [16,] "WIQDADALFGERLVAFASIEGVFFSGSFASIFWLK"              

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]          NA  0.07479708  0.05337859  0.03433479 -0.01606767 -0.04352843
       [7] -0.03729885 -0.01987896  0.04201416 -0.03060842 -0.07294796 -0.02468851
      [13]  0.01992624  0.16141861  0.02636853 -0.07032001 -0.10235492

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                                        
       [1,] "sp|P39708|DHE5_YEAST" "DTDVPAGDIGVGGR"                            
       [2,] "sp|P39708|DHE5_YEAST" "EIGYLFGAYRSYK"                             
       [3,] "sp|P07262|DHE4_YEAST" "EIGYLFGAYRSYK"                             
       [4,] "sp|P39708|DHE5_YEAST" "FHPSVNLSILK"                               
       [5,] "sp|P07262|DHE4_YEAST" "FHPSVNLSILK"                               
       [6,] "sp|P39708|DHE5_YEAST" "FHPSVNLSILKFLGFEQIFK"                      
       [7,] "sp|P07262|DHE4_YEAST" "FHPSVNLSILKFLGFEQIFK"                      
       [8,] "sp|P39708|DHE5_YEAST" "FLGFEQIFK"                                 
       [9,] "sp|P07262|DHE4_YEAST" "FLGFEQIFK"                                 
      [10,] "sp|P07262|DHE4_YEAST" "GANIASFIKVSDAMFDQGDVF"                     
      [11,] "sp|P39708|DHE5_YEAST" "GGLCVDLK"                                  
      [12,] "sp|P07262|DHE4_YEAST" "GGLCVDLK"                                  
      [13,] "sp|P39708|DHE5_YEAST" "GGLRFHPSVNLSILK"                           
      [14,] "sp|P07262|DHE4_YEAST" "GGLRFHPSVNLSILK"                           
      [15,] "sp|P39708|DHE5_YEAST" "GLNWGGSLIRPEATGFGLVYYTQAMIDYATNGKESFEGKR"  
      [16,] "sp|P39708|DHE5_YEAST" "GPYKGGLR"                                  
      [17,] "sp|P07262|DHE4_YEAST" "GPYKGGLR"                                  
      [18,] "sp|P39708|DHE5_YEAST" "IFEKFPQYK"                                 
      [19,] "sp|P39708|DHE5_YEAST" "IMINCFNDCIQAAQEYSTEK"                      
      [20,] "sp|P07262|DHE4_YEAST" "MSEPEFQQAYEEVVSSLEDSTLFEQHPEYR"            
      [21,] "sp|P39708|DHE5_YEAST" "NSWEGVLTGK"                                
      [22,] "sp|P07262|DHE4_YEAST" "NSWEGVLTGK"                                
      [23,] "sp|P39708|DHE5_YEAST" "RVTISGSGNVAQYAALK"                         
      [24,] "sp|P07262|DHE4_YEAST" "RVTISGSGNVAQYAALK"                         
      [25,] "sp|P39708|DHE5_YEAST" "SDNEIRRICYAFMR"                            
      [26,] "sp|P07262|DHE4_YEAST" "SEPEFQQAYEEVVSSLEDSTLFEQHPEYRK"            
      [27,] "sp|P39708|DHE5_YEAST" "SLEEIVDEYSTFSESK"                          
      [28,] "sp|P39708|DHE5_YEAST" "SLEEIVDEYSTFSESKMK"                        
      [29,] "sp|P07262|DHE4_YEAST" "SNNEIRRICYAFMR"                            
      [30,] "sp|P39708|DHE5_YEAST" "TSEPEFQQAYDEIVSSVEDSK"                     
      [31,] "sp|P39708|DHE5_YEAST" "TSEPEFQQAYDEIVSSVEDSKIFEKFPQYK"            
      [32,] "sp|P07262|DHE4_YEAST" "VDQELKR"                                   
      [33,] "sp|P07262|DHE4_YEAST" "VQYIAGARPWTHVQKVDIALPCATQNEVSGEEAKALVAQGVK"
      [34,] "sp|P39708|DHE5_YEAST" "VTWENDNGEQEVAQGYRVQFNSAK"                  
      [35,] "sp|P39708|DHE5_YEAST" "VTWENDNGEQEVAQGYRVQFNSAKGPYK"              
      [36,] "sp|P39708|DHE5_YEAST" "VTWTAERVDQELKK"                            

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]            NA            NA  0.0584628554 -0.1324309411 -0.0031494266
       [6] -0.0756718860  0.1882964513  0.0004911877  0.1297839721  0.1365840283
      [11] -0.0169712199 -0.0478375486  0.0897837669 -0.0939095455 -0.0177834013
      [16] -0.0239801765 -0.0587987057  0.1257093503 -0.0620754035 -0.0846196400
      [21] -0.0061624088  0.0993429132 -0.1171475105  0.0798260039  0.0339711304
      [26]  0.0372646264 -0.0126222923  0.0021501175  0.0423442955

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                    [,2]                                      
       [1,] "sp|P40212|RL13B_YEAST" "AEAEAEKKK"                               
       [2,] "sp|Q12690|RL13A_YEAST" "AEAEAEKKK"                               
       [3,] "sp|P40212|RL13B_YEAST" "AGRGFTLAEVKAAGLTAAYAR"                   
       [4,] "sp|Q12690|RL13A_YEAST" "AGRGFTLAEVKAAGLTAAYAR"                   
       [5,] "sp|P40212|RL13B_YEAST" "AISKNLPILKNHFR"                          
       [6,] "sp|Q12690|RL13A_YEAST" "AISKNLPILKNHFR"                          
       [7,] "sp|P40212|RL13B_YEAST" "APTVKYNR"                                
       [8,] "sp|Q12690|RL13A_YEAST" "APTVKYNR"                                
       [9,] "sp|P40212|RL13B_YEAST" "AREKAEAEAEK"                             
      [10,] "sp|Q12690|RL13A_YEAST" "AREKAEAEAEK"                             
      [11,] "sp|P40212|RL13B_YEAST" "AVQDNGESAFRTLR"                          
      [12,] "sp|Q12690|RL13A_YEAST" "AVQDNGESAFRTLR"                          
      [13,] "sp|P40212|RL13B_YEAST" "DGKAPEAEQVLSAAATFPIAQPATDVEAR"           
      [14,] "sp|P40212|RL13B_YEAST" "DGKAPEAEQVLSAAATFPIAQPATDVEARAVQDNGESAFR"
      [15,] "sp|P40212|RL13B_YEAST" "EKAEAEAEK"                               
      [16,] "sp|Q12690|RL13A_YEAST" "EKAEAEAEK"                               
      [17,] "sp|P40212|RL13B_YEAST" "EKAEAEAEKK"                              
      [18,] "sp|Q12690|RL13A_YEAST" "EKAEAEAEKK"                              
      [19,] "sp|Q12690|RL13A_YEAST" "EYQSKIIVFPRNGK"                          
      [20,] "sp|P40212|RL13B_YEAST" "IIVFPR"                                  
      [21,] "sp|Q12690|RL13A_YEAST" "IIVFPR"                                  
      [22,] "sp|P40212|RL13B_YEAST" "MAISKNLPILK"                             
      [23,] "sp|Q12690|RL13A_YEAST" "MAISKNLPILK"                             
      [24,] "sp|Q12690|RL13A_YEAST" "NGKAPEAEQVLSAAATFPIAQPATDVEAR"           
      [25,] "sp|P40212|RL13B_YEAST" "NLPILKNHFR"                              
      [26,] "sp|Q12690|RL13A_YEAST" "NLPILKNHFR"                              
      [27,] "sp|P40212|RL13B_YEAST" "NQEIFDANVQR"                             
      [28,] "sp|Q12690|RL13A_YEAST" "NQEIFDANVQR"                             
      [29,] "sp|P40212|RL13B_YEAST" "TIGIAVDHR"                               
      [30,] "sp|Q12690|RL13A_YEAST" "TIGIAVDHR"                               
      [31,] "sp|P40212|RL13B_YEAST" "TLRLARSEK"                               
      [32,] "sp|Q12690|RL13A_YEAST" "TLRLARSEK"                               
      [33,] "sp|P40212|RL13B_YEAST" "VSRRNAR"                                 
      [34,] "sp|Q12690|RL13A_YEAST" "VSRRNAR"                                 

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]           NA           NA -0.009031485  0.131544087 -0.034243711
       [6]  0.071698315  0.021667122  0.057868254  0.072131162 -0.050961290
      [11]  0.065263944 -0.074479885  0.017207991  0.060380335  0.106023755
      [16]  0.099660820  0.016373723  0.039791537 -0.107939659  0.021846047
      [21] -0.012865537

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                                 
       [1,] "sp|P09938|RIR2_YEAST" "AEASFWTAEEIDLSKDIHDWNNR"            
       [2,] "sp|P09938|RIR2_YEAST" "AGVMSK"                             
       [3,] "sp|P09938|RIR2_YEAST" "AYLKSHQVHRHK"                       
       [4,] "sp|P09938|RIR2_YEAST" "DIHDWNNRMNENERFFISR"                
       [5,] "sp|P09938|RIR2_YEAST" "ERTVLFPIKYHEIWQAYK"                 
       [6,] "sp|P09938|RIR2_YEAST" "LLVAFGNKK"                          
       [7,] "sp|P09938|RIR2_YEAST" "LSKDAENHKAYLK"                      
       [8,] "sp|P09938|RIR2_YEAST" "SDMLKEKLSK"                         
       [9,] "sp|P09938|RIR2_YEAST" "VENPFDFMENISLAGKTNFFEKR"            
      [10,] "sp|P09938|RIR2_YEAST" "VKSDMLKEK"                          
      [11,] "sp|P09938|RIR2_YEAST" "VSDYQK"                             
      [12,] "sp|P09938|RIR2_YEAST" "WIQDADALFGERLVAFASIEGVFFSGSFASIFWLK"
      [13,] "sp|P09938|RIR2_YEAST" "YHEIWQAYKR"                         

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]            NA  0.0439334130 -0.0679794694 -0.1319438718 -0.0375233639
       [6]  0.0677908753 -0.0163748468 -0.0538487003  0.0906672875 -0.0816404574
      [11] -0.0286334245 -0.0022983186 -0.0860850555 -0.0002294695

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
            [,1]                   [,2]                                            
       [1,] "sp|P39708|DHE5_YEAST" "ALVASGVKFVAEGANMGSTPEAISVFETARSTATNAK"         
       [2,] "sp|P39708|DHE5_YEAST" "DTDVPAGDIGVGGR"                                
       [3,] "sp|P39708|DHE5_YEAST" "ELSRHIGKDTDVPAGDIGVGGR"                        
       [4,] "sp|P07262|DHE4_YEAST" "FIAEGSNMGSTPEAIAVFETARSTATGPSEAVWYGPPK"        
       [5,] "sp|P39708|DHE5_YEAST" "FPQYKKVLPIVSVPER"                              
       [6,] "sp|P39708|DHE5_YEAST" "GGLRFHPSVNLSILK"                               
       [7,] "sp|P07262|DHE4_YEAST" "GGLRFHPSVNLSILK"                               
       [8,] "sp|P07262|DHE4_YEAST" "GLNWGGSLIRPEATGYGLVYYTQAMIDYATNGK"             
       [9,] "sp|P39708|DHE5_YEAST" "GPYKGGLR"                                      
      [10,] "sp|P07262|DHE4_YEAST" "GPYKGGLR"                                      
      [11,] "sp|P39708|DHE5_YEAST" "IFEKFPQYK"                                     
      [12,] "sp|P07262|DHE4_YEAST" "ITWTSERVDQELK"                                 
      [13,] "sp|P39708|DHE5_YEAST" "NALTGLDMGGGKGGLCVDLK"                          
      [14,] "sp|P39708|DHE5_YEAST" "RVTISGSGNVAQYAALK"                             
      [15,] "sp|P07262|DHE4_YEAST" "RVTISGSGNVAQYAALK"                             
      [16,] "sp|P39708|DHE5_YEAST" "RVTISGSGNVAQYAALKVIELGGIVVSLSDSK"              
      [17,] "sp|P39708|DHE5_YEAST" "SLEEIVDEYSTFSESK"                              
      [18,] "sp|P07262|DHE4_YEAST" "SNNEIRRICYAFMR"                                
      [19,] "sp|P39708|DHE5_YEAST" "SYKNSWEGVLTGK"                                 
      [20,] "sp|P07262|DHE4_YEAST" "SYKNSWEGVLTGK"                                 
      [21,] "sp|P07262|DHE4_YEAST" "SYKNSWEGVLTGKGLNWGGSLIRPEATGYGLVYYTQAMIDYATNGK"
      [22,] "sp|P39708|DHE5_YEAST" "TSEPEFQQAYDEIVSSVEDSK"                         
      [23,] "sp|P07262|DHE4_YEAST" "VDQELKR"                                       
      [24,] "sp|P07262|DHE4_YEAST" "VLPSLVKGANIASFIKVSDAMFDQGDVF"                  
      [25,] "sp|P07262|DHE4_YEAST" "VQYIAGARPWTHVQKVDIALPCATQNEVSGEEAKALVAQGVK"    
      [26,] "sp|P07262|DHE4_YEAST" "VQYNSAKGPYKGGLR"                               
      [27,] "sp|P39708|DHE5_YEAST" "VTWENDNGEQEVAQGYRVQFNSAK"                      
      [28,] "sp|P39708|DHE5_YEAST" "VTWTAERVDQELK"                                 

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
       [1]            NA            NA  0.0647700774  0.0569007485 -0.0384757603
       [6]  0.0095541475 -0.0063670296  0.0271016134  0.0238549146  0.0082650368
      [11]  0.0568098001 -0.0663416711 -0.0579556038 -0.1071758166 -0.0334041916
      [16] -0.0765861025 -0.0058478336  0.0204620893 -0.0578592249 -0.0379594950
      [21]  0.1119601300  0.0362848546 -0.0265005296  0.0605594207  0.0188319883
      [26] -0.0002171576

