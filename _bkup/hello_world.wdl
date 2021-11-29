task hello {
  String pattern
  File in
 
  command {
    egrep '${pattern}' '${in}'
  }
 
  runtime {
    docker: "ubuntu:16.04"
  }
 
  output {
    Array[String] matches = read_lines(stdout())
  }
}
 
workflow wf {
  call hello
}
