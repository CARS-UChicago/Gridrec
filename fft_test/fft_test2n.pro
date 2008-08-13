pro fft_test2n, data, isign
   s = size(data, /dimensions)
   nx = s[0]
   ny = s[1]
   if (!version.os eq 'linux') then object='fft_test.so' else object='Release\fft_test'
   t = call_external(object, 'fft_test2n', complex(data), nx, ny, long(isign))
end
