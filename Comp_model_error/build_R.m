function R=build_R(sig,k)
%sig is the vector cntaining 
%fourier coefficients of conductivity
%ordered as %[\sigma_0,\sigma_c,\sigm_s]=[zeroth, cosine, sine]
%k=size(NtD,1)
n=k/2;%size of NtD matrix/2
l=(length(sig)-1)/2;%length of fourier coefficients of conductivity
%[\sigma_0,\sigma_c,\sigm_s]
R_11=zeros(n,n);

for i=1:n
    for j=1:n
        if (i+j)<=l && (i+j)>0
            s=sig(i+j+1);
        else
            s=0;
        end
        if abs(i-j)==0
            ss=2*sig(1);
            m_sign=1;
        elseif abs(i-j)<=l && abs(i-j)>0
            ss=sig(abs(i-j)+1);
            m_sign=sign(i-j);
        else 
            ss=0;
            m_sign=0;
        end
            
        R_11(i,j)= s + (ss);
    end
end
R_12=zeros(n,n);

for i=1:n
    for j=1:n
        if (i+j)<=l && (i+j)>0
            s=sig(i+j+1+l);
        else
            s=0;
        end
        if abs(i-j)==0
            ss=2*sig(1);
        elseif abs(i-j)<=l && abs(i-j)>0
            ss=sig(abs(i-j)+1+l);
        else 
            ss=0;
        end
            
        R_12(i,j)=s - sign(i-j)*ss;
    end
end
R_21=zeros(n,n);

for i=1:n
    for j=1:n
        if (i+j)<=l && (i+j)>0
            s=sig(i+j+1+l);
        else
            s=0;
        end
        if abs(i-j)==0
            ss=2*sig(1);

        elseif abs(i-j)<=l && abs(i-j)>0
            ss=sig(abs(i-j)+1+l);
        else 
            ss=0;
        end
            
        R_21(i,j)= s + (sign(i-j))*ss;
    end
end
R_22=zeros(n,n);

for i=1:n
    for j=1:n
        if (i+j)<=l && (i+j)>0
            s=sig(i+j+1);
            
        else
            s=0;
        end
        if abs(i-j)==0
            ss=2*sig(1);
        elseif abs(i-j)<=l && abs(i-j)>0
            ss=sig(abs(i-j)+1);
        else 
            ss=0;
        end
            
        R_22(i,j)= -s + (ss);
    end
end
R=0.5*[R_11,R_12;R_21,R_22];
R=sparse(R);