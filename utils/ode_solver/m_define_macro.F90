#define _dealloc(a) if(allocated(a)) deallocate(a)
#define _realloc(a,n) if(allocated(a)) deallocate(a); allocate(a(n))
#define _realloc2(a,n,m) if(allocated(a)) deallocate(a); allocate(a(n,m))
#define _G write(6,*) 'GREETINGS from '//__FILE__, __LINE__
#define _T write(6,'(a,i9,3x,a18)') 'TIMINGS from '//__FILE__, __LINE__, get_cdatetime()
#define _die(m) call die(m//' in '//__FILE__, __LINE__)
#define _warn(m) call warn(m//' in '//__FILE__, __LINE__)
#define _print(m) write(6,*) m, 'm'
#define _mem_note call log_memory_note(__FILE__, iv, __LINE__)
#define _mem_note_iter(iteration) call log_memory_note(__FILE__, iv, __LINE__, iteration)
#define _zero(a) if(allocated(a)) a=0
#define _pack_indx(i,j) j*(j-1)/2+i
#define _pack_size(n) n*(n+1)/2
#define _t1 call cputime(t1)
#define _t2(t) call cputime(t2); t = t + (t2 - t1)
#define _conc(a,b) trim(a)//', '//trim(b)
