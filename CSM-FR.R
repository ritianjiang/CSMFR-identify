#The metadata of gr contains type, REF and ALT
# type represents the origin of this mutation (i.e. CEN or YC)
# REF is reference allele in hg38
# ALT is the alternative allele
gr_u<-unique(as.data.frame(sort(gr[,c("type","REF","ALT")])))

res<-data.frame()
len<-5000 #bp
i<-1
while(i<=nrow(gr_u)){
  ts<-i
  tt<-i
  sc<-c(ts)
  try({
    while(abs(gr_u[ts,]$end-gr_u[ts+1,]$end) < len){
      ts<-ts+1
      sc<-c(sc,ts)
    }
    tmp<-data.frame(seqnames=gr_u[sc[1],]$seqnames,start=gr_u[sc[1],]$start,
                    end=gr_u[sc[length(sc)],]$start,strand="*",
                    num_type=gr_u[sc,]$type %>% unique %>% sort %>% str_c(collapse = "_"),
                    locis_nums=gr_u[sc,]$start %>% unique %>% length)
  },silent = T)
  if(tmp$locis_nums>=3){res<-rbind(res,tmp)}
  i<-ts+1
}
save(res,file="Hotspots.RData")
