function result = laplacian(fs, Phi)
% LAPLACIAN Return the laplacian of a string state; [fs.i, [fs.i, Phi]]
result = commutator(fs.x, commutator(fs.x, Phi))...
    + commutator(fs.y, commutator(fs.y, Phi))...
    + commutator(fs.z, commutator(fs.z, Phi));
end