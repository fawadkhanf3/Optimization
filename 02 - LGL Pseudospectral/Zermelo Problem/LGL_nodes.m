function tau = LGL_nodes(N)

thetak = (4*[1:N]-1)*pi/(4*N+2);
sigmak = -(1-(N-1)/(8*N^3)-(39-28./sin(thetak).^2)/(384*N^4)).*cos(thetak);
    
ze = (sigmak(1:N-1)+sigmak(2:N))/2;
ep = eps*10;

ze1 = ze+ep+1;

while max(abs(ze1-ze))>=ep
    ze1 = ze;
    [dy,y] = lepoly(N,ze);
    ze = ze-(1-ze.*ze).*dy./(2*ze.*dy-N*(N+1)*y);
end

tau = [-1,ze,1]';

end