classdef BandedBatch %#ok<*PROPLC,*PREALL>

    properties ( SetAccess = private )
        n
        kl
        ku
        nbatch
        nthreads
    end

    properties ( Access = private )
        A
        ipiv
        na
        nb
    end

    methods

        function obj = BandedBatch(matrices, nthreads)

            if ( nargin == 0 || isempty(matrices) )
                return
            end

            if ( nargin == 1 )
                nthreads = maxNumCompThreads;
            end

            [m, n] = size(matrices{1});
            if (m ~= n)
                error('Matrix must be square.');
            end

            [kl, ku] = bandwidth(matrices{1});
            nbatch = length(matrices);
            lda = 2*kl+ku+1;
            na = lda*n*nbatch;
            nb =     n*nbatch;
            A    = zeros(lda, n, nbatch);
            ipiv = zeros(nb, 1);

            % Convert matrices into banded storage format for LAPACK.
            % (See https://www.netlib.org/lapack/lug/node124.html for details.)
            for k = 1:nbatch
                A(kl+1:end,:,k) = spdiags(full(matrices{k}), ku:-1:-kl).';
            end
            A = A(:);

            % Compute a banded LU decomposition of each matrix.
            mex_id_ = 'factor_banded_batch(i int, i int, i int, io double[x], o int[x], i int, i int)';
[A, ipiv] = gateway(mex_id_, n, kl, ku, A, nbatch, nthreads, na, nb);

            obj.n = n;
            obj.kl = kl;
            obj.ku = ku;
            obj.nbatch = nbatch;
            obj.nthreads = nthreads;
            obj.A = A;
            obj.ipiv = ipiv;
            obj.na = na;
            obj.nb = nb;

        end

        function x = solve(obj, b)

            n = obj.n;
            kl = obj.kl;
            ku = obj.ku;
            nbatch = obj.nbatch;
            nthreads = obj.nthreads;
            A = obj.A;
            ipiv = obj.ipiv;
            nb = obj.nb;
            na = obj.na;

            % Make a deep copy of the data so that we own the memory of x.
            x = b(:,:);

            % Solve using the banded LU decompositions.
            if ( isreal(x) )
                mex_id_ = 'solve_banded_batch(i int, i int, i int, i double[x], i int[x], io double[x], i int, i int)';
[x] = gateway(mex_id_, n, kl, ku, A, ipiv, x, nbatch, nthreads, na, nb, nb);
            else
                % If the RHS is complex then we need to do two solves:
                % one for the real part and one for the imaginary part.
                xreal = real(x);
                ximag = imag(x);
                mex_id_ = 'solve_banded_batch(i int, i int, i int, i double[x], i int[x], io double[x], i int, i int)';
[xreal] = gateway(mex_id_, n, kl, ku, A, ipiv, xreal, nbatch, nthreads, na, nb, nb);
                mex_id_ = 'solve_banded_batch(i int, i int, i int, i double[x], i int[x], io double[x], i int, i int)';
[ximag] = gateway(mex_id_, n, kl, ku, A, ipiv, ximag, nbatch, nthreads, na, nb, nb);
                x = xreal + ximag*1i;
            end

            x = reshape(x, size(b));

        end

        function x = mldivide(obj, b)
            x = solve(obj, b);
        end

    end

end
