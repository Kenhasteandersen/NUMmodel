function A3=function_convert_TM_positive(A)

N = A - diag(diag(A));
N(N>0)=0;

A2= A- N + (-N');

DA=diag(A2)'+ sum(N,1) +sum(N,2)';

A3= A2 - diag(diag(A2)) + diag(DA);

end