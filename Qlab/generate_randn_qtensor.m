function [qt] = generate_randn_qtensor(sz_tens)
%GENERATE_RANDN_QTENSOR �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    qt = quaternion(randn(sz_tens),randn(sz_tens),randn(sz_tens),randn(sz_tens));  
end

