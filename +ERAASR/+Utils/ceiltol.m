function c = ceiltol(A, tol)

    fl = floor(A);
    if A - fl < tol
        c = fl;
    else
        c = ceil(A);
    end

end