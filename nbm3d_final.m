%% set hyperparameter
h=0.1; % the unit length of a grid.
F=5;
dt=h*h/F; % CFL condition

% set a domain
x_size=2; % -x_size<=x<=x_size
y_size=2;
z_size=7;

% set the size of rectangles
x_rec = 2; % must be less than 2*x_size
y_rec = 2; % must be less than 2*y_size

% set the distance between two rectangles
dist = 5; % must be less than 2*z_size

%make a grid
x=-x_size:h:x_size;
y=-y_size:h:y_size;
z=-z_size:h:z_size;

[X,Y,Z]=meshgrid(x,y,z);

Nx=length(x);
Ny=length(y);
Nz=length(z);

phi=zeros(Nx,Ny,Nz); % level set function(before update)
phin=zeros(Nx,Ny,Nz); % level set function(after update)
zr=zeros(0,3); % zero level set of phi

%% make a sphere for making a narrow band
sq=zeros(0,3); 
for i=-2:2
    for j=-2:2
        for k=-2:2
            if(i^2+j^2+k^2<=4) % set the width of narrow band의 width by 2*h
                sq=[sq;[i,j,k]];
            end
        end
    end
end

%% Define zero level set
% make a zero level set and set function value on the outside and the inside of the zero level set
% inside = 2*h outside = -2*h

phi(:,:,:)=-2*h;
for i=1:Nx
    for j=1:Ny
        for k=1:(Nz+1)/2-dist/2/h-1
            phi(i,j,k)=2*h;
        end
        for k=(Nz+1)/2+dist/2/h+1:Nz
            phi(i,j,k)=2*h;
        end
        if(max(abs(x(i))*2/x_rec,abs(y(j))*2/y_rec)>1)
            phi(i,j,(Nz+1)/2-dist/2/h)=0;
            phi(i,j,(Nz+1)/2+dist/2/h)=0;
        elseif(max(abs(x(i))*2/x_rec,abs(y(j))*2/y_rec)==1)
            for k=(Nz+1)/2-dist/2/h:(Nz+1)/2+dist/2/h
                phi(i,j,k)=0;
            end
        else
            for k=(Nz+1)/2-dist/2/h:(Nz+1)/2+dist/2/h
                phi(i,j,k)=2*h;
            end
        end
    end
end

%% Make narrow band
[f,v]=isosurface(X,Y,Z,phi,0.0);
for n=1:length(v)
    zr=[zr;v(n,:)]; % a set of points on zero level set(but, not grid points)
end
% apporximate a point in zr into a grid point
for n=1:length(zr)
    x1=round((x_size+zr(n,1))/h+1);
    y1=round((y_size+zr(n,2))/h+1);
    z1=round((z_size+zr(n,3))/h+1);
    rx=(x_size+zr(n,1))/h+1-x1;
    ry=(y_size+zr(n,2))/h+1-y1;
    rz=(z_size+zr(n,3))/h+1-z1;
   
    % Now using a sphere with radius 2*h, set a function value on the points in narrow band such that
    % its value is the signed distance function.
    for m=1:length(sq)
        x2=sq(m,1);
        y2=sq(m,2);
        z2=sq(m,3);
        if(abs(x1+x2-round(Nx/2))<=round(Nx/2)-1 && abs(y1+y2-round(Ny/2))<=round(Ny/2)-1 && abs(z1+z2-round(Nz/2))<=round(Nz/2)-1)
            dispr=phi(x1+x2,y1+y2,z1+z2);
            dis=h*((sq(m,1)-rx)^2+(sq(m,2)-ry)^2+(sq(m,3)-rz)^2)^0.5;
            if(dispr<0)
                phi(x1+x2,y1+y2,z1+z2)=-min(dis,-dispr);
            else
                phi(x1+x2,y1+y2,z1+z2)=min(dis,dispr);
            end
        end
    end
end

%% just for convenience
ipx=zeros(1,Nx);
imx=zeros(1,Nx);
for i=1:Nx
    ipx(i)=i+1;
    imx(i)=i-1;
end
imx(1)=1;
ipx(Nx)=Nx;

ipy=zeros(1,Ny);
imy=zeros(1,Ny);
for i=1:Ny
    ipy(i)=i+1;
    imy(i)=i-1;
end
imy(1)=1;
ipy(Ny)=Ny;

ipz=zeros(1,Nz);
imz=zeros(1,Nz);
for i=1:Nz
    ipz(i)=i+1;
    imz(i)=i-1;
end
imz(1)=1;
ipz(Nz)=Nz;

%% update the value
for a=1:500 % num of iteration
    clf
    fv=isosurface(X,Y,Z,phi,0.0);
    p = patch(fv);
    set(p,'FaceColor','blue','EdgeColor','none');
    axis([-x_size x_size -y_size y_size -z_size z_size]);
    daspect([1 1 1])
    view(3);
    camlight
    lighting gouraud
    pause(.01)
    a
    
    phin=phi;
    % update the function only on the points in the narrow band
    for n=1:length(zr)
        x1=round((x_size+zr(n,1))/h+1);
        y1=round((y_size+zr(n,2))/h+1);
        z1=round((z_size+zr(n,3))/h+1);
        for m=1:length(sq)
            x2=sq(m,1);
            y2=sq(m,2);
            z2=sq(m,3);
            if(abs(x1+x2-round(Nx/2))<=round(Nx/2)-1 && abs(y1+y2-round(Ny/2))<=round(Ny/2)-1 && abs(z1+z2-round(Nz/2)<=round(Nz/2)-1))
                idx=x1+x2;
                idy=y1+y2;
                idz=z1+z2;
                pxx=(phi(ipx(idx),idy,idz)-2*phi(idx,idy,idz)+phi(imx(idx),idy,idz))/(h*h);   
                pyy=(phi(idx,ipy(idy),idz)-2*phi(idx,idy,idz)+phi(idx,imy(idy),idz))/(h*h);
                pzz=(phi(idx,idy,ipz(idz))-2*phi(idx,idy,idz)+phi(idx,idy,imz(idz)))/(h*h);
                pxy=(phi(ipx(idx),ipy(idy),idz)+phi(imx(idx),imy(idy),idz)-phi(imx(idx),ipy(idy),idz)-phi(ipx(idx),imy(idy),idz))/(4.*h*h);
                pyz=(phi(idx,ipy(idy),ipz(idz))+phi(idx,imy(idy),imz(idz))-phi(idx,imy(idy),ipz(idz))-phi(idx,ipy(idy),imz(idz)))/(4.*h*h);
                pzx=(phi(ipx(idx),idy,ipz(idz))+phi(imx(idx),idy,imz(idz))-phi(imx(idx),idy,ipz(idz))-phi(ipx(idx),idy,imz(idz)))/(4.*h*h);
                px=(phi(ipx(idx),idy,idz)-phi(imx(idx),idy,idz))/(2*h);
                py=(phi(idx,ipy(idy),idz)-phi(idx,imy(idy),idz))/(2*h);
                pz=(phi(idx,idy,ipz(idz))-phi(idx,idy,imz(idz)))/(2*h);
                curvgrad=(pxx*(py^2+px^2)+pyy*(px^2+pz^2)+pzz*(px^2+py^2)-2*pxy*px*py-2*pyz*py*pz-2*pzx*px*pz)/(2*(px^2+py^2+pz^2));
                phin(idx,idy,idz)=phi(idx,idy,idz)+curvgrad*dt;
            end
        end
    end
    
    % Now we have to make a new narrow band
    [f,v]=isosurface(X,Y,Z,phin,0.0);
    zr=zeros(0,3);
    for n=1:length(v)
        if(abs(v(n,3))<dist/2)
             zr=[zr;v(n,:)];
        end
    end
    for i=1:Nx
        for j=1:Ny
            if(max(abs(x(i))*2/x_rec,abs(y(j))*2/y_rec)>=1)
                zr=[zr;[x(i),y(j),dist/2];[x(i),y(j),-dist/2]];
                phin(i,j,(Nz+1)/2-dist/2/h)=0;
                phin(i,j,(Nz+1)/2+dist/2/h)=0;
            end
        end
    end
    
    % Edit phin
    for i=1:Nx
        for j=1:Ny
            for k=1:Nz
                if(phin(i,j,k)>0)
                    phin(i,j,k)=2*h;
                elseif(phin(i,j,k)<0)
                    phin(i,j,k)=-2*h;
                end
            end
        end
    end
    
    % Make narrow band
    for n=1:length(zr)
        x1=round((x_size+zr(n,1))/h+1);
        y1=round((y_size+zr(n,2))/h+1);
        z1=round((z_size+zr(n,3))/h+1);
        rx=(x_size+zr(n,1))/h+1-x1;
        ry=(y_size+zr(n,2))/h+1-y1;
        rz=(z_size+zr(n,3))/h+1-z1;

        for m=1:length(sq)
            x2=sq(m,1);
            y2=sq(m,2);
            z2=sq(m,3);
            if(abs(x1+x2-round(Nx/2))<=round(Nx/2)-1 && abs(y1+y2-round(Ny/2))<=round(Ny/2)-1 && abs(z1+z2-round(Nz/2))<=round(Nz/2)-1)
                dispr=phin(x1+x2,y1+y2,z1+z2);
                dis=h*((sq(m,1)-rx)^2+(sq(m,2)-ry)^2+(sq(m,3)-rz)^2)^0.5;
                if(dispr<0)
                    phin(x1+x2,y1+y2,z1+z2)=-min(dis,-dispr);
                else
                    phin(x1+x2,y1+y2,z1+z2)=min(dis,dispr);
                end
            end
        end
    end
    phi=phin;
end
