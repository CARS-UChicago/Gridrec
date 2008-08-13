pro fft_test1n, data, isign
   s = size(data, /dimensions)
   nx = s[0]
   if (!version.os eq 'linux') then object='fft_test.so' else object='Release\fft_test'
   t = call_external(object, 'fft_test1n', complex(data), nx, long(isign))
end
