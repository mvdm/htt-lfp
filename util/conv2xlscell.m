function [str] = conv2xlscell(row, cln)
    %% Convert row number & column number to work with xlswrite()
    % Example ['AAP6'] = conv2xlscell(6,718)

    
    % column
    cln2 = ceil(cln/26)-1;
    cln3 = ceil(cln2/26)-1; if cln3 < 0, cln3 = 0; end
    cln2 = cln2-cln3*26;
    cln1 = cln-cln2*26-cln3*26^2;
    vec = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

    str1 = vec(cln1);
    
    if cln2 > 0
        str2 = vec(cln2);
    else str2 = '';
    end
    
    if cln3 > 0
        str3 = vec(cln3);
    else str3 = '';
    end
                                % row
    str1 = strcat(str3,str2,str1,num2str(row));
    str = str1{1};
end

