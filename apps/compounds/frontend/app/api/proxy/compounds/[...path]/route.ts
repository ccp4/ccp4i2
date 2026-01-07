import { NextRequest, NextResponse } from 'next/server';

const DJANGO_URL = process.env.DJANGO_URL || 'http://localhost:8000';

export async function GET(
  request: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  const { path } = await params;
  const pathString = path.join('/');
  const searchParams = request.nextUrl.searchParams.toString();
  const queryString = searchParams ? `?${searchParams}` : '';
  const url = `${DJANGO_URL}/compounds/${pathString}/${queryString}`;

  console.log(`[API Proxy] GET ${url}`);

  try {
    const response = await fetch(url, {
      headers: {
        'Accept': 'application/json',
      },
    });

    const data = await response.json();
    return NextResponse.json(data, { status: response.status });
  } catch (error) {
    console.error('[API Proxy] Error:', error);
    return NextResponse.json(
      { error: 'Failed to fetch from Django' },
      { status: 500 }
    );
  }
}

export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  const { path } = await params;
  const pathString = path.join('/');
  const url = `${DJANGO_URL}/compounds/${pathString}/`;

  console.log(`[API Proxy] POST ${url}`);

  try {
    const contentType = request.headers.get('content-type') || '';
    let response: Response;

    if (contentType.includes('multipart/form-data')) {
      // Handle file uploads - pass through the FormData
      const formData = await request.formData();
      response = await fetch(url, {
        method: 'POST',
        headers: {
          'Accept': 'application/json',
          // Don't set Content-Type - let fetch set it with boundary for multipart
        },
        body: formData,
      });
    } else {
      // Handle JSON requests (or empty body for action endpoints)
      let body: any = null;
      const contentLength = request.headers.get('content-length');
      if (contentLength && parseInt(contentLength) > 0) {
        try {
          body = await request.json();
        } catch (e) {
          // Ignore JSON parse errors for empty bodies
        }
      }

      const fetchOptions: RequestInit = {
        method: 'POST',
        headers: {
          'Accept': 'application/json',
        },
      };

      if (body !== null) {
        (fetchOptions.headers as Record<string, string>)['Content-Type'] = 'application/json';
        fetchOptions.body = JSON.stringify(body);
      }

      response = await fetch(url, fetchOptions);
    }

    const data = await response.json();
    return NextResponse.json(data, { status: response.status });
  } catch (error) {
    console.error('[API Proxy] Error:', error);
    return NextResponse.json(
      { error: 'Failed to fetch from Django' },
      { status: 500 }
    );
  }
}

export async function PATCH(
  request: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  const { path } = await params;
  const pathString = path.join('/');
  const url = `${DJANGO_URL}/compounds/${pathString}/`;
  const body = await request.json();

  console.log(`[API Proxy] PATCH ${url}`);

  try {
    const response = await fetch(url, {
      method: 'PATCH',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
      },
      body: JSON.stringify(body),
    });

    const data = await response.json();
    return NextResponse.json(data, { status: response.status });
  } catch (error) {
    console.error('[API Proxy] Error:', error);
    return NextResponse.json(
      { error: 'Failed to fetch from Django' },
      { status: 500 }
    );
  }
}

export async function DELETE(
  request: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  const { path } = await params;
  const pathString = path.join('/');
  const url = `${DJANGO_URL}/compounds/${pathString}/`;

  console.log(`[API Proxy] DELETE ${url}`);

  try {
    const response = await fetch(url, {
      method: 'DELETE',
      headers: {
        'Accept': 'application/json',
      },
    });

    if (response.status === 204) {
      return new NextResponse(null, { status: 204 });
    }

    const data = await response.json();
    return NextResponse.json(data, { status: response.status });
  } catch (error) {
    console.error('[API Proxy] Error:', error);
    return NextResponse.json(
      { error: 'Failed to fetch from Django' },
      { status: 500 }
    );
  }
}
