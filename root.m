
fun = @f;
z = fzero(fun,0)


function y=f(flux)
    %out = 17/8*flux - 15/8*abs(flux);
    y = 30*flux - flux.^2/2 + sign(flux-20).*(flux-20).^2/2 - sign(flux-40).*(flux-40).^2/2 - 25*flux;

end
