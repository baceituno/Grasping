function gray = grayCodes(base,nbit)
% Gray Code Generation

seqlength = base^nbit;
code = zeros(seqlength,nbit);

for i = 1:nbit
    RepetaionRange= base^(nbit-i+1)/base;
    IncrementSeq = 0;
    sequenceElement = 1;
    while(IncrementSeq < seqlength)
        end_point = IncrementSeq+RepetaionRange;
        start_point = IncrementSeq+1;
        code(start_point:end_point,i) = sequenceElement - 1;
        IncrementSeq = IncrementSeq + RepetaionRange;
        sequenceElement = sequenceElement + 1;
        if(sequenceElement > base)
            sequenceElement = 1;
        end
    end
end

gray = zeros(seqlength,nbit);

for j = 1:seqlength
    g(1) = code(j,1);
    for i = 2:nbit
        x = xor(code(j,i-1), code(j,i));
        g(i) = x;
    end
    gray(j,:) = g(:);
end