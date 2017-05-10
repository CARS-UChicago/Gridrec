pro fft_test2f, data, isign
   common fft_test_common, library

   s = size(data, /dimensions)
   nx = s[0]
   ny = s[1]
   t = call_external(library, 'fft_test2f', complex(data), nx, ny, long(isign))
end
