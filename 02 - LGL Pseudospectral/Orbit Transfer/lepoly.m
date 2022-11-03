function [varargout] = lepoly(n,x)

if nargout == 1
   if n == 0, varargout{1} = ones(size(x)); return; end
   if n == 1, varargout{1} = x; return; end
   polylst = ones(size(x)); 
   poly = x;
   
   for k = 2:n
       polyn = ((2*k-1)*x.*poly-(k-1)*polylst)/k;
       polylst = poly;
       poly = polyn;
   end
   varargout{1} = polyn;
end
   
if nargout == 2
    
    if n == 0, varargout{2} = ones(size(x)); varargout{1} = zeros(size(x)); return; end
    if n == 1, varargout{2} = x; varargout{1} = ones(size(x)); return; end
    
    polylst = ones(size(x));
    pderlst = zeros(size(x));
    poly = x;
    pder = ones(size(x));
    
    for k = 2:n
        polyn = ((2*k-1)*x.*poly-(k-1)*polylst)/k;
        pdern = pderlst+(2*k-1)*poly;
        
        polylst = poly;
        poly = polyn;
        
        pderlst = pder;
        pder = pdern;
    end
    varargout{2} = polyn;
    varargout{1} = pdern;
   
end