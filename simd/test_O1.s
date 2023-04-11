	.file	"test.cpp"
	.text
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC0:
	.string	"%f\n"
.LC1:
	.string	"%d\n"
	.text
	.globl	main
	.type	main, @function
main:
.LFB6930:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$8, %rsp
	.cfi_offset 13, -24
	.cfi_offset 12, -32
	.cfi_offset 3, -40
	movl	$4000000000, %eax
	subq	%rax, %rsp
	movq	%rsp, %r13
	subq	%rax, %rsp
	movq	%rsp, %r12
	subq	%rax, %rsp
	leaq	3(%rsp), %rbx
	movq	%rbx, %rdx
	shrq	$2, %rdx
	andq	$-4, %rbx
	movl	$0, %eax
.L2:
	movss	0(%r13,%rax,4), %xmm0
	addss	(%r12,%rax,4), %xmm0
	movss	%xmm0, (%rbx,%rax,4)
	addq	$1, %rax
	cmpq	$1000000000, %rax
	jne	.L2
	pxor	%xmm0, %xmm0
	cvtss2sd	0(,%rdx,4), %xmm0
	leaq	.LC0(%rip), %rdi
	movl	$1, %eax
	call	printf@PLT
	movq	%r13, %rsi
	andl	$63, %esi
	leaq	.LC1(%rip), %rdi
	movl	$0, %eax
	call	printf@PLT
	movq	%r12, %rsi
	andl	$63, %esi
	leaq	.LC1(%rip), %rdi
	movl	$0, %eax
	call	printf@PLT
	movq	%rbx, %rsi
	andl	$63, %esi
	leaq	.LC1(%rip), %rdi
	movl	$0, %eax
	call	printf@PLT
	movl	$0, %eax
	leaq	-24(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE6930:
	.size	main, .-main
	.type	_GLOBAL__sub_I_main, @function
_GLOBAL__sub_I_main:
.LFB7423:
	.cfi_startproc
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	leaq	_ZStL8__ioinit(%rip), %rdi
	call	_ZNSt8ios_base4InitC1Ev@PLT
	leaq	__dso_handle(%rip), %rdx
	leaq	_ZStL8__ioinit(%rip), %rsi
	movq	_ZNSt8ios_base4InitD1Ev@GOTPCREL(%rip), %rdi
	call	__cxa_atexit@PLT
	addq	$8, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE7423:
	.size	_GLOBAL__sub_I_main, .-_GLOBAL__sub_I_main
	.section	.init_array,"aw"
	.align 8
	.quad	_GLOBAL__sub_I_main
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
	.hidden	__dso_handle
	.ident	"GCC: (Debian 10.2.1-6) 10.2.1 20210110"
	.section	.note.GNU-stack,"",@progbits
