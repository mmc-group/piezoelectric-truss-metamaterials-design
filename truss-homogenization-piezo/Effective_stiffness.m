function [C,e,mu]=Effective_stiffness(T,KK,Sdof,xx,yy,zz,Fnodes_full)
% % Load case:1
X011=zeros(Sdof,1);
X011(1:7:end)=xx;
f11=KK*X011;% % Equivalent Nodal Forces
ff11=T'*f11;
KKK=T'*KK*T;
Xg11=pinv(KKK)*ff11;
X1_11=T*Xg11;
FF_11=KK*X1_11;
C11=((X011-X1_11)'*(f11-FF_11));

% % Load case:2
X022=zeros(Sdof,1);
X022(2:7:end)=yy;
f22=KK*X022;    % % Equivalent Nodal Forces
ff22=T'*f22;
Xg22=pinv(KKK)*ff22;
X1_22=T*Xg22;
FF_22=KK*X1_22;
C22=((X022-X1_22)'*(f22-FF_22));

% % Load case:3
X033=zeros(Sdof,1);
X033(3:7:end)=zz;
f33=KK*X033;    % % Equivalent Nodal Forces
ff33=T'*f33;
Xg33=pinv(KKK)*ff33;
X1_33=T*Xg33;
FF_33=KK*X1_33;
C33=((X033-X1_33)'*(f33-FF_33));

% % Load case:4
X044=zeros(Sdof,1);
X044(3:7:end)=0.5*yy;
X044(2:7:end)=0.5*zz;
f44=KK*X044;    % % Equivalent Nodal Forces
ff44=T'*f44;
Xg44=pinv(KKK)*ff44;
X1_44=T*Xg44;
FF_44=KK*X1_44;
C44=((X044-X1_44)'*(f44-FF_44));% C32=C23;

% % Load case:5
X055=zeros(Sdof,1);
X055(3:7:end)=0.5*xx;
X055(1:7:end)=0.5*zz;
f55=KK*X055;    % % Equivalent Nodal Forces
ff55=T'*f55;
Xg55=pinv(KKK)*ff55;
X1_55=T*Xg55;
FF_55=KK*X1_55;
C55=((X055-X1_55)'*(f55-FF_55));% C31=C13;

% % Load case:6
X066=zeros(Sdof,1);
X066(2:7:end)=0.5*xx;
X066(1:7:end)=0.5*yy;
f66=KK*X066;    % % Equivalent Nodal Forces
ff66=T'*f66;
Xg66=pinv(KKK)*ff66;
X1_66=T*Xg66;
FF_66=KK*X1_66;
C66=((X066-X1_66)'*(f66-FF_66));% C21=C12;

% % Load case 7
X077=zeros(Sdof,1);
X077(7:7:end)=-xx;
f77=KK*X077;    % % Equivalent Nodal Forces
ff77=T'*f77;
Xg77=pinv(KKK)*ff77;
X1_77=T*Xg77;
FF_77=KK*X1_77;
mu11=((X077-X1_77)'*(f77-FF_77));% mu11=C77;

% % Load case 8
X088=zeros(Sdof,1);
X088(7:7:end)=-yy;
f88=KK*X088;    % % Equivalent Nodal Forces
ff88=T'*f88;
Xg88=pinv(KKK)*ff88;
X1_88=T*Xg88;
FF_88=KK*X1_88;
mu22=((X088-X1_88)'*(f88-FF_88));% mu22=C88;

% % Load case 9
X099=zeros(Sdof,1);
X099(7:7:end)=-zz;
f99=KK*X099;    % % Equivalent Nodal Forces
ff99=T'*f99;
Xg99=pinv(KKK)*ff99;
X1_99=T*Xg99;
FF_99=KK*X1_99;
mu33=((X099-X1_99)'*(f99-FF_99));% mu33=C99;

% % Off-diagonal terms

C12=((X011-X1_11)'*(f22-FF_22));
C13=((X011-X1_11)'*(f33-FF_33));
C14=((X011-X1_11)'*(f44-FF_44));
C15=((X011-X1_11)'*(f55-FF_55));
C16=((X011-X1_11)'*(f66-FF_66));
e11=((X077-X1_77)'*(f11-FF_11)); %e11=C71
e12=((X077-X1_77)'*(f22-FF_22)); %e12=C72
e13=((X077-X1_77)'*(f33-FF_33)); %e13=C73
e14=((X077-X1_77)'*(f44-FF_44)); %e14=C74
e15=((X077-X1_77)'*(f55-FF_55)); %e15=C75
e16=((X077-X1_77)'*(f66-FF_66)); %e16=C76
e21=((X088-X1_88)'*(f11-FF_11)); %e21=C81
e22=((X088-X1_88)'*(f22-FF_22)); %e22=C82
e23=((X088-X1_88)'*(f33-FF_33)); %e23=C83
e24=((X088-X1_88)'*(f44-FF_44)); %e24=C84
e25=((X088-X1_88)'*(f55-FF_55)); %e25=C85
e26=((X088-X1_88)'*(f66-FF_66)); %e25=C86
e31=((X099-X1_99)'*(f11-FF_11)); %e31=C91
e32=((X099-X1_99)'*(f22-FF_22)); %e32=C92
e33=((X099-X1_99)'*(f33-FF_33)); %e33=C93
e34=((X099-X1_99)'*(f44-FF_44)); %e34=C94
e35=((X099-X1_99)'*(f55-FF_55)); %e35=C95
e36=((X099-X1_99)'*(f66-FF_66)); %e36=C96
mu12=((X077-X1_77)'*(f88-FF_88));% mu33=C99;
mu13=((X077-X1_77)'*(f99-FF_99));% mu33=C99;
mu21=((X088-X1_88)'*(f77-FF_77));% mu33=C99;
mu23=((X088-X1_88)'*(f99-FF_99));% mu33=C99;
mu31=((X099-X1_99)'*(f77-FF_77));% mu33=C99;
mu32=((X088-X1_88)'*(f99-FF_99));% mu33=C99;



C21=((X022-X1_22)'*(f11-FF_11)); 
C23=((X022-X1_22)'*(f33-FF_33));
C24=((X022-X1_22)'*(f44-FF_44));
C25=((X022-X1_22)'*(f55-FF_55));
C26=((X022-X1_22)'*(f66-FF_66));
C31=((X033-X1_33)'*(f11-FF_11));
C32=((X033-X1_33)'*(f22-FF_22));
C34=((X033-X1_33)'*(f44-FF_44));
C35=((X033-X1_33)'*(f55-FF_55));
C36=((X033-X1_33)'*(f66-FF_66));
C41=((X044-X1_44)'*(f11-FF_11));
C42=((X044-X1_44)'*(f22-FF_22));
C43=((X044-X1_44)'*(f33-FF_33));
C45=((X044-X1_44)'*(f55-FF_55));
C46=((X044-X1_44)'*(f66-FF_66));
C51=((X055-X1_55)'*(f11-FF_11));
C52=((X055-X1_55)'*(f22-FF_22));
C53=((X055-X1_55)'*(f33-FF_33));
C54=((X055-X1_55)'*(f44-FF_44));
C56=((X055-X1_55)'*(f66-FF_66));
C61=((X066-X1_66)'*(f11-FF_11));
C62=((X066-X1_66)'*(f22-FF_22));
C63=((X066-X1_66)'*(f33-FF_33));
C64=((X066-X1_66)'*(f44-FF_44));
C65=((X066-X1_66)'*(f55-FF_55));
C=[C11 C12 C13 C14 C15 C16
   C21 C22 C23 C24 C25 C26
   C31 C32 C33 C34 C35 C36
   C41 C42 C43 C44 C45 C46
   C51 C52 C53 C54 C55 C56
   C61 C62 C63 C64 C65 C66];

e=-[e11 e12  e13 e14 e15 e16
   e21 e22 e23 e24 e25 e26
   e31 e32 e33 e34 e35 e36];

mu=[mu11 mu12 mu13
    mu21 mu22 mu23
    mu31 mu32 mu33];

end