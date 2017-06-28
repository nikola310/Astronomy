format long;
%konstante
Kr = 1.67e-23;     % specificna refraktivnost CO2 u cm3/molec
rv = 6052;        % poluprecnik venere
ra = 6220;        % visina do vrha atmosfere
ds = 108208000;   % udaljenost sunce-venera
dsz = 1.496e8;      % sunce-zemlja
s = 1.392684e6;     % poluprecnik sunca
zvs = pi - asin(dsz*sin(atan(s/dsz))/ds); % ugao zemlja-venera-sunce
dz = dsz*sin(pi-zvs-atan(s/dsz))/sin(zvs);  % udaljenost venera-zemlja
%   rezolucija
pixel = 0.5;
sloj = 0.5;
% granice x i y
x_1 = -6200;   x_2 = -6100;
y_1 = pixel;   y_2 = 2;

for x_osa = x_1:pixel:x_2;
    for y_osa = y_1:pixel:y_2;
        if sqrt(x_osa^2 + y_osa^2) > 6125 && sqrt(x_osa^2 + y_osa^2) < 6170
            % racunanje centralnog zraka
            % tacka na oreolu
            x = x_osa + pixel/2;
            y = y_osa + pixel/2;
            % uglovi pod kojima zrak sa zemlje prodire u venerinu atmosferu
            alpha0 = atan(sqrt(x^2+y^2)/dz);
            fi = acos(x/sqrt(x^2+y^2));
            % deo pre petlje se odnosi na tacku prodora u prvom sloju atmosfere
            % upadni i izlazni ugao
            n1 = 1;
            n2 = 1 + Kr*10^(-0.1*(ra-rv)+25);
            alpha1 = asin(dz*sin(alpha0)/ra); %upadni ugao
            beta1 = asin(n1*sin(alpha1)/n2);          %izlazni ugao
            % koordinate tacke prodora kroz sloj u ravni prelamanja zraka
            % gde je x osa paralelna sa zrakom sa zemlje, a (0,0) je u centru venere
            y1 = ra*sin(alpha1);
            x1 = -ra*cos(alpha1);
            % ugao pod kojim je prelomljeni zrak odstupio od upadnog zraka
            % i racunanje parametara njegove jednacine prave
            delta = alpha1 - beta1;
            k = tan(pi-delta);
            c = y1 - k*x1;
            % provera da li prelomljeni zrak ima presek sa narednim slojem
            r = ra - sloj;
            D = (k^2*r^2)+r^2-c^2; %ako je diskriminanta D>=0 onda postoji presek
            %b=0;
            while D >= 0 % petlja se odnosi na tacku prodora za trenutni sloj
                %b=b+1
                % racunanje koordinate prodora kroz sloj
                x2 = (-k*c - sqrt(D))/(k^2+1);
                y2 = k*x2+c;
                plot(x2,y2); axis equal;
                hold on;
                % racunanje upadnog i izlaznog ugla iz sloja
                n1 = n2;   %indeks prelamanja prethodnog sloja
                n2 = 1 + Kr*10^(-0.1*(r-rv)+25);   %indeks prelamanja narednog sloja
                gamma2 = acos(x2/r);
                alpha2 = pi - gamma2 - delta;
                beta2 = asin(n1*sin(alpha2)/n2);
                % ugao pod kojim je prelomljeni zrak odstupio od upadnog zraka
                % i racunanje parametara njegove jednacine prave
                delta2 = alpha2 - beta2;
                delta = delta + delta2;
                k = tan(pi-delta);
                c = y2 - k*x2;
                % provera da li prelomljeni zrak ima presek sa narednim slojem
                r = r - sloj;
                D = (k^2*r^2)+r^2-c^2; %ako je diskriminanta D>=0 onda postoji presek
            end
            % racunanje parametara tacke izlaska zraka (Ta) iz venerine atmosfere
            c3 = y1 - k*x1;
            x3 = (-k*c3 + sqrt(k^2*ra^2+ra^2-c3^2))/(k^2+1);
            y3 = k*x3+c3;
            gamma = asin(y3/ra);
            % pomocni uglovi izmedju razlicitih pravaca u prostoru
            gamma3 = gamma + alpha0;                                              % ugao polozaja Ta u odnosu na pravac zemlja-venera
            delta0 = acos(-cos(zvs)*cos(gamma3) - sin(zvs)*sin(gamma3)*cos(fi));  % ugao Ta-venera-sunce
            mi = atan(ra*sin(delta0)/(ds-ra*cos(delta0)));                        % ugao venera-sunce-Ta
            tau = acos((-cos(zvs)-cos(delta0)*cos(gamma3))/(sin(delta0)*sin(gamma3)));
            if gamma3 >= alpha1
                epsilon = acos(cos(delta0)*cos(alpha1) + sin(delta0)*sin(alpha1)*cos(tau));
                epsilonk = acos(cos(delta0+mi)*cos(alpha1) + sin(delta0+mi)*sin(alpha1)*cos(tau));
                ksi1 = acos((cos(alpha1)-cos(epsilon)*cos(delta0))/(sin(epsilon)*sin(delta0)));
                ksi2 = acos((cos(gamma3)+cos(delta0)*cos(zvs))/(sin(delta0)*sin(zvs))) - ksi1;
                rs = sqrt(ra^2*sin(delta0)^2 + ((ds-ra*cos(delta0))*tan(epsilon))^2 + 2*ra*sin(delta0)*(ds-ra*cos(delta0))*tan(epsilon)*cos(ksi1));
                fi_s = ksi2 + asin(ra*sin(delta0)*sin(ksi1)/rs);
            else
                epsilon = acos(-cos(alpha1-gamma3)*cos(zvs)+sin(alpha1-gamma3)*sin(zvs)*cos(fi));
                ksi1 = acos((cos(gamma3)+cos(zvs)*cos(delta0))/(sin(delta0)*sin(zvs)));
                ksi2 = acos((cos(alpha1-gamma3)+cos(epsilon)*cos(zvs))/(sin(epsilon)*sin(zvs)));
                kappa = pi - ksi1 - ksi2;
                epsilonk = acos(cos(mi)*cos(epsilon)+sin(mi)*sin(epsilon)*cos(kappa));
                rs = sqrt(ra^2*sin(delta0)^2 + ((ds-ra*cos(delta0))*tan(epsilon))^2 + 2*ra*sin(delta0)*(ds-ra*cos(delta0))*tan(epsilon)*cos(ksi1+ksi2));
                fi_s = ksi1 - asin((ds-ra*cos(delta0))*tan(epsilon)*sin(ksi1+ksi2)/rs);
            end
            xs0 = rs*cos(fi_s);
            ys0 = rs*sin(fi_s);
            tas = sqrt(ra^2 + ds^2 - 2*ra*ds*cos(delta0));
            if epsilonk <= atan(s/tas)   % ako centralni zrak pogadja sunce
                % racunanje 4 tacke
                for x_pix = 0:1;
                    for y_pix = 0:1;
                        % tacka na oreolu
                        x = x_osa + pixel*x_pix;
                        y = y_osa + pixel*y_pix;
                        % uglovi pod kojima zrak sa zemlje prodire u venerinu atmosferu
                        alpha0 = atan(sqrt(x^2+y^2)/dz);
                        fi = acos(x/sqrt(x^2+y^2));
                        % deo pre petlje se odnosi na tacku prodora u prvom sloju atmosfere
                        % upadni i izlazni ugao
                        n1 = 1;
                        n2 = 1 + Kr*10^(-0.1*(ra-rv)+25);
                        alpha1 = asin(dz*sin(alpha0)/ra);         %upadni ugao
                        beta1 = asin(n1*sin(alpha1)/n2);          %izlazni ugao
                        % koordinate tacke prodora kroz sloj u ravni prelamanja zraka
                        % gde je x osa paralelna sa zrakom sa zemlje, a (0,0) je u centru venere
                        y1 = ra*sin(alpha1);
                        x1 = -ra*cos(alpha1);
                        % ugao pod kojim je prelomljeni zrak odstupio od upadnog zraka
                        % i racunanje parametara njegove jednacine prave
                        delta = alpha1 - beta1;
                        k = tan(pi-delta);
                        c = y1 - k*x1;
                        % provera da li prelomljeni zrak ima presek sa narednim slojem
                        r = ra - sloj;
                        D = (k^2*r^2)+r^2-c^2; %ako je diskriminanta D>=0 onda postoji presek
                        while D >= 0 % petlja se odnosi na tacku prodora za trenutni sloj
                            % racunanje koordinate prodora kroz sloj
                            x2 = (-k*c - sqrt(D))/(k^2+1);
                            y2 = k*x2+c;
                            %plot(x2,y2); axis equal;
                            %hold on;
                            % racunanje upadnog i izlaznog ugla iz sloja
                            n1 = n2;   %indeks prelamanja prethodnog sloja
                            n2 = 1 + Kr*10^(-0.1*(r-rv)+25);   %indeks prelamanja narednog sloja
                            gamma2 = acos(x2/r);
                            alpha2 = pi - gamma2 - delta;
                            beta2 = asin(n1*sin(alpha2)/n2);
                            % ugao pod kojim je prelomljeni zrak odstupio od upadnog zraka
                            % i racunanje parametara njegove jednacine prave
                            delta2 = alpha2 - beta2;
                            delta = delta + delta2;
                            k = tan(pi-delta);
                            c = y2 - k*x2;
                            % provera da li prelomljeni zrak ima presek sa narednim slojem
                            r = r - sloj;
                            D = (k^2*r^2)+r^2-c^2; %ako je diskriminanta D>=0 onda postoji presek
                        end
                        % racunanje parametara tacke izlaska zraka (Ta) iz venerine atmosfere
                        c3 = y1 - k*x1;
                        x3 = (-k*c3 + sqrt(k^2*ra^2+ra^2-c3^2))/(k^2+1);
                        y3 = k*x3+c3;
                        gamma = asin(y3/ra);
                        % pomocni uglovi izmedju razlicitih pravaca u prostoru
                        gamma3 = gamma + alpha0;                                              % ugao polozaja Ta u odnosu na pravac zemlja-venera
                        delta0 = acos(-cos(zvs)*cos(gamma3) - sin(zvs)*sin(gamma3)*cos(fi));  % ugao Ta-venera-sunce
                        mi = atan(ra*sin(delta0)/(ds-ra*cos(delta0)));                        % ugao venera-sunce-Ta
                        tau = acos((-cos(zvs)-cos(delta0)*cos(gamma3))/(sin(delta0)*sin(gamma3)));
                        if gamma3 >= alpha1
                            epsilon = acos(cos(delta0)*cos(alpha1) + sin(delta0)*sin(alpha1)*cos(tau));
                            epsilonk = acos(cos(delta0+mi)*cos(alpha1) + sin(delta0+mi)*sin(alpha1)*cos(tau));
                            ksi1 = acos((cos(alpha1)-cos(epsilon)*cos(delta0))/(sin(epsilon)*sin(delta0)));
                            ksi2 = acos((cos(gamma3)+cos(delta0)*cos(zvs))/(sin(delta0)*sin(zvs))) - ksi1;
                            rs = sqrt(ra^2*sin(delta0)^2 + ((ds-ra*cos(delta0))*tan(epsilon))^2 + 2*ra*sin(delta0)*(ds-ra*cos(delta0))*tan(epsilon)*cos(ksi1));
                            fi_s = ksi2 + asin(ra*sin(delta0)*sin(ksi1)/rs);
                        else
                            epsilon = acos(-cos(alpha1-gamma3)*cos(zvs)+sin(alpha1-gamma3)*sin(zvs)*cos(fi));
                            ksi1 = acos((cos(gamma3)+cos(zvs)*cos(delta0))/(sin(delta0)*sin(zvs)));
                            ksi2 = acos((cos(alpha1-gamma3)+cos(epsilon)*cos(zvs))/(sin(epsilon)*sin(zvs)));
                            kappa = pi - ksi1 - ksi2;
                            epsilonk = acos(cos(mi)*cos(epsilon)+sin(mi)*sin(epsilon)*cos(kappa));
                            rs = sqrt(ra^2*sin(delta0)^2 + ((ds-ra*cos(delta0))*tan(epsilon))^2 + 2*ra*sin(delta0)*(ds-ra*cos(delta0))*tan(epsilon)*cos(ksi1+ksi2));
                            fi_s = ksi1 - asin((ds-ra*cos(delta0))*tan(epsilon)*sin(ksi1+ksi2)/rs);
                        end
                        tas = sqrt(ra^2 + ds^2 - 2*ra*ds*cos(delta0));  
                        if epsilonk > atan(s/tas)
                            xs = s*cos(fi_s);
                            ys = s*sin(fi_s);
                        else
                            xs = rs*cos(fi_s);
                            ys = rs*sin(fi_s);
                        end
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
                %                 diag_1 = sqrt((xs1-xs0)^2+(ys1-ys0)^2);
                %                 diag_2 = sqrt((xs2-xs0)^2+(ys2-ys0)^2);
                %                 diag_3 = sqrt((xs3-xs0)^2+(ys3-ys0)^2);
                %                 diag_4 = sqrt((xs4-xs0)^2+(ys4-ys0)^2);
                str_a = sqrt((xs1-xs2)^2+(ys1-ys2)^2);
                str_b = sqrt((xs2-xs3)^2+(ys2-ys3)^2);
                str_c = sqrt((xs3-xs4)^2+(ys3-ys4)^2);
                str_d = sqrt((xs4-xs1)^2+(ys4-ys1)^2);
                diag = sqrt((xs1-xs3)^2+(ys1-ys3)^2);
                sp1 = (str_a+str_b+diag)/2;
                sp2 = (diag+str_c+str_d)/2;
                triangle1 = sqrt(sp1*(sp1-str_a)*(sp1-str_b)*(sp1-diag));
                triangle2 = sqrt(sp2*(sp2-diag)*(sp2-str_c)*(sp2-str_d));
                Sv = triangle1 + triangle2;
                sigma = asin(sin(epsilonk)*tas/s) - epsilonk;
                zrak = s*sin(sigma)/sin(epsilonk);
                z = sqrt(zrak^2+ra^2+2*zrak*ra*cos(alpha1));
                sigma0 = atan(rs/(ds-z*cos(asin(rs/z))));
                Sn = Sv/cos(sigma0);
                Sk = Sn*cos(sigma);      
                Slika((x_osa-x_1)/pixel+1, y_osa/pixel+1) = Sk;
            else
                Slika((x_osa-x_1)/pixel+1, y_osa/pixel+1) = 0;
            end
        else
            Slika((x_osa-x_1)/pixel+1, y_osa/pixel+1) = 0;
        end
    end
end