function [ B, C ] = SKPDecomposition( A, SizeB, SizeC, Hermitian )
%SKPDECOMPOSITION 
% sum of Kronecker product decomposition
    m = size( A, 1 );
    n = size( A, 2 );
    m1 = SizeB( 1 );
    n1 = SizeB( 2 );
    m2 = SizeC( 1 );
    n2 = SizeC( 2 );
    
    if nargin < 4
        Hermitian = false;
    end

    assert( m == m1 * m2 );
    assert( n == n1 * n2 );
    
    if Hermitian
        assert( m1 == n1 );
        assert( m2 == n2 );
        A = 0.5 * ( A + A' );
    end
    
    R = reshape( permute( reshape( A, [ m2, m1, n2, n1 ] ), [ 2 4 1 3 ] ), m1 * n1, m2 * n2 );
    [ U, S, V ] = svd(R);
    M = m1*n1;
    if m1*n1 > m2*n2
        M = m2*n2;
    end
    B=cell(1,M);
    C=cell(1,M);
    for ii=1:M
        SqrtS = sqrt( S(ii,ii) );
        B1 = reshape( U(:,ii) * SqrtS, m1, n1 );
        C1 = reshape( V(:,ii) * SqrtS, m2, n2 );
        B1 = B1*real(C1(1,1));
        C1 = C1/real(C1(1,1));
        if Hermitian
            B1 = 0.5 * ( B1 + B1' );
            C1 = 0.5 * ( C1 + C1' );
            if all( diag( B1 ) < 0 ) && all( diag( C1 ) < 0 )
                B1 = -B1;
                C1 = -C1;
            end
        end
        B{1,ii}=B1;
        C{1,ii}=conj(C1);
    end
end

