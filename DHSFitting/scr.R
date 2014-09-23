library('Rcpp')


cppFunction('IntegerVector examp2(Rcpp::List activeList) {
int n = activeList.size();
IntegerVector out(10);
IntegerVector numactive(n);
   for(int tt = 0; tt < n; ++tt) {
        IntegerVector active = activeList[tt];
for(int ii = 0; ii < active.size(); ++ii) {
        Rprintf(" %d", active(ii));
        }
    Rprintf("\\n");
    for(IntegerVector::iterator jj = active.begin(); jj != active.end(); ++jj) { 
    Rprintf(" %d", *jj);
        }
    Rprintf("\\n");
}
  return out;
}')

numeric_limits<int>::min()

eList = list(1:10, c(1:6,8), 1:6, 2:3, 3, 1)
alf <- examp2(eList)

cppFunction('double ex3(double x) {
  return -INFINITY;
}')


cppFunction('bool ex3(LogicalVector x) {
bool out = is_true(any(x));
  return out;
}')

ex3(c(F,F,F))

cppFunction('bool anyC(LogicalVector x) {
  return x.any();
}')


anyC(c(2,3))


cppFunction('double ex2(double x) {
double z=3/x;
if (x == std::isnan()) {
Rprintf("err");
z = 0;
}
  return z;
}')

ex2(2)
ex2(0)



args(precLoop)
