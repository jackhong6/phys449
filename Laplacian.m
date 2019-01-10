function result = laplacian(fs, Phi)
% LAPLACIAN Return the laplacian of a string state
    N = size(fs.x, 1);
    
    A = - (N^2 - 1) / (4 * fs.R^4);

    result = A * (commutator(fs.x, commutator(fs.x, Phi))...
        + commutator(fs.y, commutator(fs.y, Phi))...
        + commutator(fs.z, commutator(fs.z, Phi)));
end