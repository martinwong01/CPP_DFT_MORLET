	.file	"testsimd.cpp"
	.text
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC0:
	.string	"%f\n"
	.section	.text.startup,"ax",@progbits
	.p2align 4
	.globl	main
	.type	main, @function
main:
.LFB5368:
	.cfi_startproc
	leaq	8(%rsp), %r10
	.cfi_def_cfa 10, 0
	andq	$-32, %rsp
	movl	$4000000000, %eax
	pushq	-8(%r10)
	pushq	%rbp
	.cfi_escape 0x10,0x6,0x2,0x76,0
	movq	%rsp, %rbp
	pushq	%r10
	.cfi_escape 0xf,0x3,0x76,0x78,0x6
	subq	$8, %rsp
	subq	%rax, %rsp
	movq	%rsp, %rsi
	subq	%rax, %rsp
	movq	%rsp, %rcx
	subq	%rax, %rsp
	xorl	%eax, %eax
	leaq	3(%rsp), %rdx
	movq	%rdx, %rdi
	andq	$-4, %rdx
	shrq	$2, %rdi
	.p2align 4,,10
	.p2align 3
.L2:
	vmovups	(%rsi,%rax,4), %xmm2
	vmovups	(%rcx,%rax,4), %xmm3
	vinsertf128	$0x1, 16(%rsi,%rax,4), %ymm2, %ymm0
	vinsertf128	$0x1, 16(%rcx,%rax,4), %ymm3, %ymm1
	vaddps	%ymm1, %ymm0, %ymm0
	vmovups	%xmm0, (%rdx,%rax,4)
	vextractf128	$0x1, %ymm0, 16(%rdx,%rax,4)
	addq	$8, %rax
	cmpq	$1000000000, %rax
	jne	.L2
	vxorps	%xmm0, %xmm0, %xmm0
	movl	$1, %eax
	vcvtss2sd	0(,%rdi,4), %xmm0, %xmm0
	leaq	.LC0(%rip), %rdi
	vzeroupper
	call	printf@PLT
	movq	-8(%rbp), %r10
	.cfi_def_cfa 10, 0
	xorl	%eax, %eax
	leave
	leaq	-8(%r10), %rsp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5368:
	.size	main, .-main
	.ident	"GCC: (Debian 10.2.1-6) 10.2.1 20210110"
	.section	.note.GNU-stack,"",@progbits
