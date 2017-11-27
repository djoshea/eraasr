function f = floortol(A, tol)

    cl = ceil(A);
    if cl - A < tol
        f = cl;
    else
        f = floor(A);
    end

end