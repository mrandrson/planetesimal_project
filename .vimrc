set autoindent
set smartindent
set tabstop=3
set shiftwidth=3
set ruler
set incsearch
if has("gui_running")
   if has("gui_gtk2")
      set guifont=Bitstream\ Vera\ Sans\ Mono\ 14
   elseif has("x11")
      set guifont=-*-Bitstream\ Vera\ Sans\ Mono-medium-r-normal-*-*-180-*-*-m-*-*
   else
      set guifont=Bitstream\ Vera\ Sans\ Mono:h14:cDEFAULT
   endif
endif 
set backup
