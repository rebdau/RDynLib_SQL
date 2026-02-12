Overall.net_SQL <- function(sql_path,
                            exp.id,
                            dbkey,
                            nr_of_seq = 2,
                            thr1 = 3,
                            thr2 = 0.1,
                            thr3 = 0.2) {
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  par(mfrow = c(1,2), mar = c(1,1,1,1))
  
  ##CSPP
  net.lst1 <- net.addit_SQL(sql_path = sql_path,
                            exp.id = exp.id,
                            nettype = "CSPP",
                            thr1 = thr1,
                            thr2 = thr2,
                            thr3 = thr3)
  
  locNET1.list <- Local.net_SQL(net.lst1,
                                sql_path = sql_path,
                                dbkey = dbkey,
                                nettype = "CSPP",
                                nr_of_seq = nr_of_seq)
  
  
  ##GNPS 
  net.lst2 <- net.addit_SQL(sql_path = sql_path,
                            exp.id = exp.id,
                            nettype = "GNPS",
                            thr1 = thr1,
                            thr2 = thr2,
                            thr3 = thr3)
  
  locNET2.list <- Local.net_SQL(net.lst2,
                                sql_path = sql_path,
                                dbkey = dbkey,
                                nettype = "GNPS",
                                nr_of_seq = nr_of_seq)
  
  return(list(locNET1.list, locNET2.list))
}
