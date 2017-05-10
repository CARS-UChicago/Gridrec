pro fft_test1i, data, isign
   common fft_test_common, library

   s = size(data, /dimensions)
   nx = s[0]
   t = call_external(library, 'fft_test1i', complex(data), nx, long(isign))
end
