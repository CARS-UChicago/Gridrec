pro fft_test1d, nx=nx, isign=isign, nloop=nloop, f0, f1, f2, f3   
   common fft_test_common, library

   if (!version.os eq 'linux') then library='O.linux-x86_64/libfftTestIDL.so' else library='O.windows-x64/fftTestIDL.dll'

   if (n_elements(nx) eq 0)    then nx = 2048
   if (n_elements(nloop) eq 0) then nloop = 2048
   if (n_elements(isign) eq 0) then isign=1
   data = complexarr(nx)
   data[200:300]=complex(1.5,.5)
   t0 = systime(1)
   for i=0, nloop-1 do begin
      f0 = fft(data, isign)
   endfor
   ; Renormalize
   if (isign eq -1) then f0 = f0 * nx
   print, 'Time with IDL internal FFT', systime(1)-t0

   t0 = systime(1)
   for i=0, nloop-1 do begin
      f1 = data
      fft_test1n, f1, isign
   endfor
   print, 'Time with Numerical Recipes FFT', systime(1)-t0

   t0 = systime(1)
   for i=0, nloop-1 do begin
      f3 = data
      fft_test1f, f3, isign
   endfor
   print, 'Time with FFTW', systime(1)-t0

   print, 'Differences:'
   diff10 = f1 - f0;
   print, ' Numerical Recipes and IDL:  max= ', $
           max(abs(diff10)), ' RMS = ', sqrt(total(abs(diff10)^2)/nx)
   diff30 = f3 - f0;
   print, ' FFTW and IDL:               max= ', $
           max(abs(diff30)), ' RMS = ', sqrt(total(abs(diff30)^2)/nx)
   diff31 = f3 - f1;
   print, ' FFTW and Numerical Recipes: max= ', $
           max(abs(diff31)), ' RMS = ', sqrt(total(abs(diff31)^2)/nx)
end
