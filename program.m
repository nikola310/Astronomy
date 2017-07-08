clc, clear all;
format long;
%figure;

%konstante
Kr = 1.671145065169295e-23;     % specificna refraktivnost CO2 u cm3/molec
rv = 6052;        % poluprecnik venere
ra = 6252;        % visina do vrha atmosfere
ds = 108208000;   % udaljenost sunce-venera
dsz = 1.496e8;      % sunce-zemlja
s = 1.392684e6;     % poluprecnik sunca
zvs = pi - asin(dsz*sin(atan(s/dsz))/ds); % ugao zemlja-venera-sunce                    0.9832*
dz = dsz*sin(pi-zvs-atan(s/dsz))/sin(zvs);  % udaljenost venera-zemlja                  0.9832*

% parametri funkcije za sin_gamma
a_0 =   0.968010519226915;      %2.98981e-5;
b0 =   1.59947101441224e-4;
A1 =   -0.00136023580793088;
xc =   60;
t1 =   1.15122600695638;
A2 =   -0.00201713322711304;
t2 =   3.09415632368965;
A3 =   -0.00224254394960332;
t3 =   4.56319313649013;

%parametri funkcije za Rejlijevo rasejanje
a0  =  -11.94264894519979;
a1  =  0.02514911522051399;
a2  =  0.00429223317264581;
a3  =  2.679999606151944E-5;
a4  =  -1.219219409753835E-7;
a5  =  -7.389905968244192E-9;
a6  =  -7.826908176623142E-11;
b1  =  0.0504392303349202;
b2  =  -0.00378817105831951;
b3  =  4.090459666544743E-5;
b4  =  4.442767437496151E-7;
b5  =  -4.464139326408204E-9;
b6  =  1.837412375053375E-11;

%   rezolucija
pixel = 0.5;

% granice x i y
x_1 = 6052;   x_2 = 6252;
y_1 = pixel;   y_2 = 12;

a=0;  %brojac petlje

Slika = zeros((x_2-x_1)/pixel+1, (y_2-y_1)/pixel+1);

for x_osa = x_1:pixel:x_2;
%     a = a+1;
%     if a/1000 == round(a/1000)
%         a
%     end

    %b=0;
        
    for y_osa = y_1:pixel:y_2;
        %         b=b+1
        %         if b/10 == round(b/10)
        %             b
        %         end
        % petlja za odgovarajuci prsten
        
        
        if sqrt(x_osa^2 + y_osa^2) > 6122 && sqrt(x_osa^2 + y_osa^2) < 6160
            
            promasaj = 0;
            
            % racunanje centralnog zraka
            
            % tacka na oreolu
            x = x_osa + pixel/2;
            y = y_osa + pixel/2;
            
            % uglovi pod kojima zrak sa zemlje prodire u venerinu atmosferu
            r0 = sqrt(x^2+y^2);
            alpha0 = atan(r0/dz);
            fi = acos(x/r0);
            alpha1 = asin(dz*sin(alpha0)/ra);  alpha1_0 = alpha1; %upadni ugao
            
            % funkcija
            % nezavisna promenljiva
            r = r0 - 6052; %r sloja - r venere
            
            % parametri funkcije na pocetku
            
            sin_gamma = a_0+b0*r + A1*exp(-(r-xc)/t1)+A2*exp(-(r-xc)/t2)+A3*exp(-(r-xc)/t3);
             
            gamma = asin(sin_gamma);  
            
            % pomocni uglovi izmedju razlicitih pravaca u prostoru
            gamma3 = gamma + alpha0;                                              % ugao polozaja Ta u odnosu na pravac zemlja-venera
            delta0 = acos(-cos(zvs)*cos(gamma3) - sin(zvs)*sin(gamma3)*cos(fi));  % ugao Ta-venera-sunce
            mi = atan(ra*sin(delta0)/(ds-ra*cos(delta0)));                        % ugao venera-sunce-Ta
            tau = acos((-cos(zvs)-cos(delta0)*cos(gamma3))/(sin(delta0)*sin(gamma3)));
            epsilon = acos(cos(delta0)*cos(alpha1) + sin(delta0)*sin(alpha1)*cos(tau)); 
            kat1 = ra*sin(delta0);
            kat2 = (ds-ra*cos(delta0))*tan(epsilon);
                
            if gamma3 >= alpha1
                epsilonk0 = acos(cos(delta0+mi)*cos(alpha1) + sin(delta0+mi)*sin(alpha1)*cos(tau));
                ksi1 = acos((cos(alpha1)-cos(epsilon)*cos(delta0))/(sin(epsilon)*sin(delta0))); 
                ksi2 = acos((cos(gamma3)+cos(delta0)*cos(zvs))/(sin(delta0)*sin(zvs))) - ksi1;
                rs = sqrt(kat1^2 + kat2^2 + 2*kat1*kat2*cos(ksi1));  rs0=rs;
                if kat1^2 <= kat2^2 + rs^2
                    fi_s = ksi2 + asin(kat1*sin(ksi1)/rs);            %     fi_s0 = fi_s;
                else
                    fi_s = ksi2 + pi - asin(kat1*sin(ksi1)/rs);       %     fi_s0 = fi_s;
                end
            else
                epsilonk0 = acos(cos(delta0+mi)*cos(alpha1) + sin(delta0+mi)*sin(alpha1)*cos(tau));
                ksi1 = acos((cos(gamma3)+cos(zvs)*cos(delta0))/(sin(delta0)*sin(zvs)));
                ksi2 = acos((cos(alpha1)-cos(delta0)*cos(epsilon))/(sin(delta0)*sin(epsilon))) - ksi1; 
                rs = sqrt(kat1^2 + kat2^2 + 2*kat1*kat2*cos(ksi1+ksi2));  rs0=rs;
                if kat2^2 <= kat1^2 + rs^2
                    fi_s = ksi1 - asin(kat2*sin(ksi1+ksi2)/rs);          %   fi_s0 = fi_s;
                else
                    fi_s = ksi1 - pi + asin(kat2*sin(ksi1+ksi2)/rs);     %   fi_s0 = fi_s;
                end
            end
            
            
            xs0 = rs*cos(fi_s);
            ys0 = rs*sin(fi_s);
            
            tas0 = sqrt(ra^2 + ds^2 - 2*ra*ds*cos(delta0));
            
            if epsilonk0 <= tan(s/tas0)   % ako centralni zrak pogadja sunce
                
                
                % racunanje 4 tacke
                
                for x_pix = 0:1;
                    for y_pix = 0:1;
                        % tacka na oreolu
                        x = x_osa + pixel*x_pix;
                        y = y_osa + pixel*y_pix;
                        
                        
                        % uglovi pod kojima zrak sa zemlje prodire u venerinu atmosferu
                        r0 = sqrt(x^2+y^2);
                        alpha0 = atan(r0/dz);
                        fi = acos(x/r0);
                        alpha1 = asin(dz*sin(alpha0)/ra);   %upadni ugao
                        
                        % funkcija
                        % nezavisna promenljiva
                        r = r0 - 6052;
                        
                        % parametri funkcije na pocetku
                        
                        sin_gamma = a_0+b0*r+ A1*exp(-(r-xc)/t1)+A2*exp(-(r-xc)/t2)+A3*exp(-(r-xc)/t3);
                        
                        gamma = asin(sin_gamma);
                        
                        % pomocni uglovi izmedju razlicitih pravaca u prostoru
                        gamma3 = gamma + alpha0;                                              % ugao polozaja Ta u odnosu na pravac zemlja-venera
                        delta0 = acos(-cos(zvs)*cos(gamma3) - sin(zvs)*sin(gamma3)*cos(fi));  % ugao Ta-venera-sunce
                        mi = atan(ra*sin(delta0)/(ds-ra*cos(delta0)));                        % ugao venera-sunce-Ta
                        tau = acos((-cos(zvs)-cos(delta0)*cos(gamma3))/(sin(delta0)*sin(gamma3)));
                        epsilon = acos(cos(delta0)*cos(alpha1) + sin(delta0)*sin(alpha1)*cos(tau));
                        kat1 = ra*sin(delta0);
                        kat2 = (ds-ra*cos(delta0))*tan(epsilon);
                            
                        if gamma3 >= alpha1
                            epsilonk = acos(cos(delta0+mi)*cos(alpha1) + sin(delta0+mi)*sin(alpha1)*cos(tau));
                            ksi1 = acos((cos(alpha1)-cos(epsilon)*cos(delta0))/(sin(epsilon)*sin(delta0)));
                            ksi2 = acos((cos(gamma3)+cos(delta0)*cos(zvs))/(sin(delta0)*sin(zvs))) - ksi1;
                            rs = sqrt(kat1^2 + kat2^2 + 2*kat1*kat2*cos(ksi1));
                            if kat1^2 <= kat2^2 + rs^2
                                fi_s = ksi2 + asin(kat1*sin(ksi1)/rs);
                            else
                                fi_s = ksi2 + pi - asin(kat1*sin(ksi1)/rs);
                            end
                        else
                            epsilonk = acos(cos(delta0+mi)*cos(alpha1) + sin(delta0+mi)*sin(alpha1)*cos(tau));
                            ksi1 = acos((cos(gamma3)+cos(zvs)*cos(delta0))/(sin(delta0)*sin(zvs)));
                            ksi2 = acos((cos(alpha1)-cos(delta0)*cos(epsilon))/(sin(delta0)*sin(epsilon))) - ksi1;
                            rs = sqrt(kat1^2 + kat2^2 + 2*kat1*kat2*cos(ksi1+ksi2));
                            if kat2^2 <= kat1^2 + rs^2
                                fi_s = ksi1 - asin(kat2*sin(ksi1+ksi2)/rs); 
                            else
                                fi_s = ksi1 - pi + asin(kat2*sin(ksi1+ksi2)/rs);
                            end
                        end
                        
                        tas = sqrt(ra^2 + ds^2 - 2*ra*ds*cos(delta0));
                        
                        if epsilonk > atan(s/tas)
                            promasaj = promasaj + 1;
                        end
                        
                        xs = rs*cos(fi_s);
                        ys = rs*sin(fi_s);
                        
                        
                        if x_pix == 0 && y_pix == 0
                            xs1 = xs;
                            ys1 = ys;
                        elseif x_pix == 1 && y_pix == 0
                            xs2 = xs;
                            ys2 = ys;
                        elseif x_pix == 1 && y_pix == 1
                            xs3 = xs;
                            ys3 = ys;
                        elseif x_pix == 0 && y_pix == 1
                            xs4 = xs;
                            ys4 = ys;
                        end
                                                
                    end
                end
                
%                 plot(xs1, ys1, 'o', xs2, ys2, 'o', xs3, ys3, 'o', xs4, ys4, 'o')
%                 hold on;
                
                
                %         diag_1 = sqrt((xs1-xs0)^2+(ys1-ys0)^2);
                %         diag_2 = sqrt((xs2-xs0)^2+(ys2-ys0)^2);
                %         diag_3 = sqrt((xs3-xs0)^2+(ys3-ys0)^2);
                %         diag_4 = sqrt((xs4-xs0)^2+(ys4-ys0)^2);
                
                str_a = sqrt((xs1-xs2)^2+(ys1-ys2)^2);
                str_b = sqrt((xs2-xs3)^2+(ys2-ys3)^2);
                str_c = sqrt((xs3-xs4)^2+(ys3-ys4)^2);
                str_d = sqrt((xs4-xs1)^2+(ys4-ys1)^2);
                diag  = sqrt((xs1-xs3)^2+(ys1-ys3)^2);
                
                sp1 = (str_a+str_b+diag)/2;
                sp2 = (diag+str_c+str_d)/2;
                
                triangle1 = sqrt(sp1*(sp1-str_a)*(sp1-str_b)*(sp1-diag)); 
                triangle2 = sqrt(sp2*(sp2-diag)*(sp2-str_c)*(sp2-str_d));
                
                Sv = triangle1 + triangle2;
                
                % korekcije sjaja:
                % 1. Limb darkening
                
                sigma_Ta = asin(sin(epsilonk0)*tas0/s) - epsilonk0;   % sa tacke Ta
                mi_Ta = cos(sigma_Ta);
                
                zrak = s*sin(sigma_Ta)/sin(epsilonk0);
                z = sqrt(zrak^2+ra^2+2*zrak*ra*cos(alpha1_0));
                
                sigma_V = asin(z*sin(atan(rs0/ds))/s);    % sa Venere
                mi_V = cos(sigma_V);
                
                Sn = Sv/(1-0.587*(1-mi_V^0.5)+0.199*(1-mi_V)-0.823*(1-mi_V^1.5)+0.433*(1-mi_V^2));
                Sk = Sn*(1-0.587*(1-mi_Ta^0.5)+0.199*(1-mi_Ta)-0.823*(1-mi_Ta^1.5)+0.433*(1-mi_Ta^2))*(5-promasaj)/5;
                
                % 2. Rejlijevo rasejanje
                % parametri na pocetku
                
                log_apsorp = (a0 + a1*r + a2*r^2 + a3*r^3 + a4*r^4 + a5*r^5 + a6*r^6)/(1 + b1*r + b2*r^2 + b3*r^3 + b4*r^4 + b5*r^5 + b6*r^6);
                trans = 1 - 10^log_apsorp;   % transparencija
                
                Sk = Sk * trans;  %konacni intenzitet pravougaonika
                
                Slika(round((x_osa-x_1)/pixel)+1, round(y_osa/pixel)+1) = Sk;
           
            else
                Slika(round((x_osa-x_1)/pixel)+1, round(y_osa/pixel)+1) = 0;
            end
                
        else
            Slika(round((x_osa-x_1)/pixel)+1, round(y_osa/pixel)+1) = 0;
        end
        
        
    end
    
end
 
save('Venera2.mat', 'Slika')
Slika1 = (real(Slika'));   %.^(1/2.3);
% amin = min(Slika1(:));
% amax = max(Slika1(:));
% Slika1 = (Slika1/amax);
imagesc(Slika1); colormap(gray); axis equal;