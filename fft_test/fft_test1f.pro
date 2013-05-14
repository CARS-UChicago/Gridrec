pro fft_test1f, data, isign
   common fft_test_common, library

   s = size(data, /dimensions)
   nx = s[0]
   t = call_external(library, 'fft_test1f', complex(data), nx, long(isign))
end
