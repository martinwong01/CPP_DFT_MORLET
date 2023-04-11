	.file	"test.cpp"
	.text
	.section	.rodata
	.type	_ZStL19piecewise_construct, @object
	.size	_ZStL19piecewise_construct, 1
_ZStL19piecewise_construct:
	.zero	1
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
.LC0:
	.string	"%f\n"
.LC1:
	.string	"%d\n"
	.text
	.globl	main
	.type	main, @function
main:
.LFB5662:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$88, %rsp
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	movq	%rsp, %rax
	movq	%rax, %rbx
	movq	$1000000000, -72(%rbp)
	movq	-72(%rbp), %rax
	leaq	-1(%rax), %rdx
	movq	%rdx, -64(%rbp)
	movq	%rdx, %rax
	addq	$1, %rax
	movq	%rax, -128(%rbp)
	movq	$0, -120(%rbp)
	movq	%rdx, %rax
	addq	$1, %rax
	movq	%rax, %r14
	movl	$0, %r15d
	movq	%rdx, %rax
	addq	$1, %rax
	leaq	0(,%rax,4), %rdx
	movl	$16, %eax
	subq	$1, %rax
	addq	%rdx, %rax
	movl	$16, %ecx
	movl	$0, %edx
	divq	%rcx
	imulq	$16, %rax, %rax
	subq	%rax, %rsp
	movq	%rsp, %rax
	addq	$3, %rax
	shrq	$2, %rax
	salq	$2, %rax
	movq	%rax, -56(%rbp)
	movq	-72(%rbp), %rax
	subq	$1, %rax
	movq	%rax, -88(%rbp)
	movq	%rax, %rdx
	addq	$1, %rdx
	movq	%rdx, %r12
	movl	$0, %r13d
	movq	%rax, %rdx
	addq	$1, %rdx
	movq	%rdx, %r10
	movl	$0, %r11d
	addq	$1, %rax
	leaq	0(,%rax,4), %rdx
	movl	$16, %eax
	subq	$1, %rax
	addq	%rdx, %rax
	movl	$16, %ecx
	movl	$0, %edx
	divq	%rcx
	imulq	$16, %rax, %rax
	subq	%rax, %rsp
	movq	%rsp, %rax
	addq	$3, %rax
	shrq	$2, %rax
	salq	$2, %rax
	movq	%rax, -96(%rbp)
	movq	-72(%rbp), %rax
	subq	$1, %rax
	movq	%rax, -104(%rbp)
	movq	%rax, %rdx
	addq	$1, %rdx
	movq	%rdx, %r8
	movl	$0, %r9d
	movq	%rax, %rdx
	addq	$1, %rdx
	movq	%rdx, %rsi
	movl	$0, %edi
	addq	$1, %rax
	leaq	0(,%rax,4), %rdx
	movl	$16, %eax
	subq	$1, %rax
	addq	%rdx, %rax
	movl	$16, %ecx
	movl	$0, %edx
	divq	%rcx
	imulq	$16, %rax, %rax
	subq	%rax, %rsp
	movq	%rsp, %rax
	addq	$3, %rax
	shrq	$2, %rax
	salq	$2, %rax
	movq	%rax, -112(%rbp)
	movq	$0, -80(%rbp)
	movq	$0, -80(%rbp)
.L3:
	movq	-80(%rbp), %rax
	cmpq	-72(%rbp), %rax
	jnb	.L2
	movq	-56(%rbp), %rax
	movq	-80(%rbp), %rdx
	movss	(%rax,%rdx,4), %xmm1
	movq	-96(%rbp), %rax
	movq	-80(%rbp), %rdx
	movss	(%rax,%rdx,4), %xmm0
	addss	%xmm1, %xmm0
	movq	-112(%rbp), %rax
	movq	-80(%rbp), %rdx
	movss	%xmm0, (%rax,%rdx,4)
	addq	$1, -80(%rbp)
	jmp	.L3
.L2:
	movq	-112(%rbp), %rax
	movss	(%rax), %xmm0
	pxor	%xmm2, %xmm2
	cvtss2sd	%xmm0, %xmm2
	movq	%xmm2, %rax
	movq	%rax, %xmm0
	leaq	.LC0(%rip), %rdi
	movl	$1, %eax
	call	printf@PLT
	movq	-56(%rbp), %rax
	andl	$63, %eax
	movq	%rax, -72(%rbp)
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	leaq	.LC1(%rip), %rdi
	movl	$0, %eax
	call	printf@PLT
	movq	-96(%rbp), %rax
	andl	$63, %eax
	movq	%rax, -72(%rbp)
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	leaq	.LC1(%rip), %rdi
	movl	$0, %eax
	call	printf@PLT
	movq	-112(%rbp), %rax
	andl	$63, %eax
	movq	%rax, -72(%rbp)
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	leaq	.LC1(%rip), %rdi
	movl	$0, %eax
	call	printf@PLT
	movq	%rbx, %rsp
	movl	$0, %eax
	leaq	-40(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5662:
	.size	main, .-main
	.type	_Z41__static_initialization_and_destruction_0ii, @function
_Z41__static_initialization_and_destruction_0ii:
.LFB6155:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movl	%edi, -4(%rbp)
	movl	%esi, -8(%rbp)
	cmpl	$1, -4(%rbp)
	jne	.L7
	cmpl	$65535, -8(%rbp)
	jne	.L7
	leaq	_ZStL8__ioinit(%rip), %rdi
	call	_ZNSt8ios_base4InitC1Ev@PLT
	leaq	__dso_handle(%rip), %rdx
	leaq	_ZStL8__ioinit(%rip), %rsi
	movq	_ZNSt8ios_base4InitD1Ev@GOTPCREL(%rip), %rax
	movq	%rax, %rdi
	call	__cxa_atexit@PLT
.L7:
	nop
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE6155:
	.size	_Z41__static_initialization_and_destruction_0ii, .-_Z41__static_initialization_and_destruction_0ii
	.type	_GLOBAL__sub_I_main, @function
_GLOBAL__sub_I_main:
.LFB6156:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movl	$65535, %esi
	movl	$1, %edi
	call	_Z41__static_initialization_and_destruction_0ii
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE6156:
	.size	_GLOBAL__sub_I_main, .-_GLOBAL__sub_I_main
	.section	.init_array,"aw"
	.align 8
	.quad	_GLOBAL__sub_I_main
	.hidden	__dso_handle
	.ident	"GCC: (Debian 10.2.1-6) 10.2.1 20210110"
	.section	.note.GNU-stack,"",@progbits
