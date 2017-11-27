function tf = isequaltol(A, B, tol)

tf = abs(A-B) < tol;

end