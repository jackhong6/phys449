function result = laplacian(J, Phi)
    result = commutator(J.x, commutator(J.x, Phi))...
           + commutator(J.y, commutator(J.y, Phi))...
           + commutator(J.z, commutator(J.z, Phi));
end

