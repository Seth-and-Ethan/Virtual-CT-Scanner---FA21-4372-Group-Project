function Image = SheppLogan(phantom_img, img_size, NT, Anginc, Probe)

theta = 1:Anginc:180; 
PhantomPixelSize = 0.3;

% Transducer = Linear
if Probe == "Linear"
    Length = 20; % The length of the probe
    if NT > 1
        xpitch = round(Length/(NT-1)*(1/PhantomPixelSize)); 
        ypitch = 0;
    else
        xpitch = 0; ypitch = 0;
    end
% Transducer = Arc
elseif Probe == "Arc"
    Radius = 200; % The radius of the probe in mm
    Arc = 45; % The total angle of probe in degrees
    if NT > 1 
        angpitch = Arc/(NT-1);
        xpitch = round(2*cosd(angpitch/2)*Radius*(1/PhantomPixelSize)); 
        ypitch = round(2*sind(angpitch/2)*Radius*(1/PhantomPixelSize)); 
    else
       angpitch = 0; xpitch = 0; ypitch = 0;
    end 
end
if xpitch > ypitch
    maxpitch = xpitch;
else
    maxpitch = ypitch;
end
if NT > 1
    maxoffset = (round((maxpitch *(NT-1)))); % The maximum offset of the phantom in pixels
else
    maxoffset = 0;
end

PhantomSize = img_size; % The dimensions of the phantom in pixels 
org_phantom = phantom_img; % Modified Shepp-Logan phantom
gen = zeros(PhantomSize+maxoffset); % Create PhantomSize+maxoffset x PhantomSize+maxoffset matrix 

% Reconstructs the phantom based on the probe chosen by the user
Image = zeros(img_size); 
if Probe == "Linear"
    gen(round(ypitch+1):round(ypitch)+PhantomSize,round(xpitch+1):round(xpitch)+img_size) = org_phantom;
    for t=1:NT
        xStart = round((t)*xpitch*(PhantomPixelSize));
        yStart = round((t)*ypitch*(PhantomPixelSize));
        pTemp = gen(yStart+1:yStart+PhantomSize,xStart+1:xStart+PhantomSize);
        [r,~] = radon(pTemp,theta);
        t1=iradon(r,theta);
        mid = round(length(t1)/2);
        Image = Image + t1(mid-(img_size/2)+1:mid+(img_size/2),mid-(img_size/2)+1:mid+(img_size/2));
    end
elseif Probe == "Arc"
    Angle = angpitch;
    gen = org_phantom;
    for t=1:NT
        pTemp = imrotate(gen,(t-1)*Angle,'nearest');
        [r,~] = radon(pTemp,theta);
        t1 = iradon(r,theta);
        t1 = imrotate(t1,-(t-1)*Angle,'nearest');
        mid = round(length(t1)/2);
        Image = Image + t1(mid-(img_size/2)+1:mid+(img_size/2),mid-(img_size/2)+1:mid+(img_size/2));
    end
end
Image = Image/NT; % return the reconstructed image / number of transducers, as this will make the image clearer
end