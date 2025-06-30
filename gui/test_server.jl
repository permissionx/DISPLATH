using HTTP

# Include logging utilities
include("../src/logging.jl")

function simple_handler(req::HTTP.Request)
    if req.target == "/"
        return HTTP.Response(200, ["Content-Type" => "text/html"], """
        <html>
        <head><title>DISPLAΘ Test</title></head>
        <body>
            <h1>🚀 DISPLAΘ BCA Simulator</h1>
            <p>Test server is working!</p>
            <p>Time: $(now())</p>
        </body>
        </html>
        """)
    else
        return HTTP.Response(404, ["Content-Type" => "text/plain"], "Not found")
    end
end

log_success("Starting test server on port 8081...")
HTTP.serve(simple_handler, "0.0.0.0", 8081)